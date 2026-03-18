from django.http import JsonResponse
import os
import pandas as pd
import json
import time
from ggawes import settings
import datetime
import subprocess
from datetime import datetime, timedelta
import uuid
from django.db import connection
import io
from decimal import Decimal
import psutil
import subprocess
from zipfile import ZipFile
import glob
from django.contrib.auth.hashers import make_password, check_password
from django.core.mail import send_mail
from Crypto.Cipher import AES
import base64
from collections import Counter
import sqlite3

#########################################setting####################################################
NS2000_path = '/bioinfo/NS2000/RD/' #NS2000_path chipid位置
nextflow="/opt/nextflow"
# fail_mail_list=['ChenGyiLin@BionetTX.com','KennethYang@GGA.ASIA']
# success_mail_list=['ChenGyiLin@BionetTX.com','KennethYang@GGA.ASIA']
fail_mail_list=['ChenGyiLin@BionetTX.com']
success_mail_list=['ChenGyiLin@BionetTX.com']
KEY = 'GGANIPTSYSTEM0123456789876543210'.encode()
#########################################function###################################################
def sqlexe(query, params=None):
    try:
        with connection.cursor() as cursor:
            cursor.execute(query, params)
            connection.commit()
    except Exception as e:
        print("Database error:", e)

def sqlquery(query, params=None):
    try:
        with connection.cursor() as cursor:
            cursor.execute(query, params)
            columns = [col[0] for col in cursor.description]
            rows = cursor.fetchall()
            return pd.DataFrame(rows, columns=columns)
    except Exception as e:
        pass
        return pd.DataFrame()

def convertTime(df, column_name):
    df[column_name] = pd.to_datetime(df[column_name], errors='coerce')
    df[column_name] = df[column_name].apply(lambda x: x.strftime('%Y/%m/%d') if pd.notnull(x) else '-')
    return df

def coerce(v):
    if isinstance(v, Decimal):
        f = float(v)
        return int(f) if f.is_integer() else f
    return v

def run_qc_to_database(result_path):
    runqc_file = os.path.join(result_path, "run_qc.csv")
    runqc_df = pd.read_csv(runqc_file, encoding="utf-8", dtype=str)
    rename_map = {
            'chipid': 'chipSNo',
            'Average_%>=Q30': 'average_q30',
            'Yield_Gbp': 'yield_gbp',
            '%Clusters_PF': 'clusters_ph',
            '%Occupied': 'occupied'
        }
    runqc_df = runqc_df.rename(columns=rename_map)
    row = runqc_df.iloc[0]
    chip_sno        = row.get('chipSNo')
    avg_q30_pct     = row.get('average_q30')
    clusters_pf_pct = row.get('clusters_ph')
    occupied_pct    = row.get('occupied')
    yield_gbp       = row.get('yield_gbp')

    sqlexe(
        """
        INSERT INTO run_qc ("chipSNo", "average_q30", "clusters_ph", "occupied", "yield_gbp")
        VALUES (%s, %s, %s, %s, %s)
        ON CONFLICT ("chipSNo") DO NOTHING
        """,
        [chip_sno, avg_q30_pct, clusters_pf_pct, occupied_pct, yield_gbp]
    )

def sample_qc_to_database(result_path):
    sampleqc_file = os.path.join(result_path, "sample_qc.csv")
    sampleqc_df = pd.read_csv(sampleqc_file, encoding="utf-8", dtype=str)

    rename_map = {
        "chipid": "chipSNo",
        "sample": "sampleSNo",
        "sex": "sex",
        "input_reads": "input_reads",
        "%_>=10x": "10x",
        "mean_target_coverage": "mean_target_coverage",
        "uniformity_of_coverage": "uniformity_of_coverage",
        "aligned_reads": "aligned_reads",
        "%_aligned": "aligned",
        "%_enrichment": "enrichment",
        "%_padded_enrichment": "padded_enrichment",
    }
    sampleqc_df = sampleqc_df.rename(columns=rename_map)
    sampleqc_df["UniID"] = sampleqc_df["sampleSNo"].astype(str).str.strip() + "_" + sampleqc_df["chipSNo"].astype(str).str.strip()

    sql = f"""INSERT INTO "QCanalysis" ("UniID","sampleSNo","chipSNo","sex",
            "input_reads","10x","mean_target_coverage","uniformity_of_coverage",
            "aligned_reads","aligned","enrichment","padded_enrichment") 
            VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)
            ON CONFLICT ("UniID") DO NOTHING
        """

    cols_order = [
        "UniID","sampleSNo","chipSNo","sex",
        "input_reads","10x","mean_target_coverage","uniformity_of_coverage",
        "aligned_reads","aligned","enrichment","padded_enrichment"
    ]
    sampleqc_df = sampleqc_df[cols_order]
    for index, row in sampleqc_df.iterrows():
        data = tuple(map(str, row))
        sqlexe(sql, data)
        
def process_result_summary(result_path, chipid):
    sql = f"""UPDATE "sample_info" SET "status" = %s WHERE "chipSNo" = %s"""
    sqlexe(sql, ["分析完成", chipid])
    samplelist = os.listdir(result_path)
    run_qc_to_database(result_path)
    sample_qc_to_database(result_path)
    variant_to_database(result_path,chipid)

def check_run_completion(folder_path: str) -> bool:
    """
    先檢查 CopyComplete.txt；若無，再檢查 Analysis/1/Data/report.html。
    任一存在即視為完成，回傳 True；否則 False。
    """
    # 1) CopyComplete.txt
    copy_complete_path = os.path.join(folder_path, "CopyComplete.txt")
    if os.path.isfile(copy_complete_path):
        return True

    # 2) Analysis/1/Data/report.html（次要判斷）
    report_path = os.path.join(folder_path, "Analysis", "1", "Data", "report.html")
    if os.path.isfile(report_path):
        return True

    return False

def update_chip_info(chipid, status):
    chipsql = f"""UPDATE "chip_info" SET "status" = %s WHERE "chipSNo" = %s"""
    sqlexe(chipsql, [status, chipid])

def sending_mail_fail(chipid,output_path):
    subject = '[通知]WES Nextflow 分析失敗'
    message = f'Chip ID: {chipid} 的 Nextflow 分析失敗，請查看{output_path}檢查相關問題。'
    recipient_list = fail_mail_list
    send_mail(subject, message, settings.EMAIL_HOST_USER, recipient_list)

def sending_mail_success(chipid):
    subject = '[通知]WES Nextflow 分析完成'
    message = f'Chip ID: {chipid} 的 Nextflow 分析完成，請由上機管理進行查詢'
    recipient_list = success_mail_list
    send_mail(subject, message, settings.EMAIL_HOST_USER, recipient_list)

def import_to_database(output_path):
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(f"{output_path}/database.txt", "w") as f:
        f.write(f"Import to database completed on {current_time}.\n")

def get_sequencing_info(chipid):
    if not os.path.exists(NS2000_path):
        return None, [], 0
    folder_list = [d for d in os.listdir(NS2000_path) if os.path.isdir(os.path.join(NS2000_path, d))]
    for folder_name in folder_list:
        parts = folder_name.split("_")
        if len(parts) < 2:
            continue
        # 取日期
        sequencing_date = None
        try:
            sequencing_date = datetime.strptime(parts[0], "%y%m%d").strftime("%Y-%m-%d")
        except ValueError:
            pass
        folder_chipid = parts[-1]
        if folder_chipid != chipid:
            continue
        folder_path = os.path.join(NS2000_path, folder_name)
        
        # 先確認 CopyComplete.txt
        if not check_run_completion(folder_path):
            print(f"{chipid} sequencing not finish!!")
            return sequencing_date, [], 0

        # 讀 SampleSheet
        sample_sheet_path = os.path.join(folder_path, "SampleSheet.csv")
        sample_ids = []
        if os.path.exists(sample_sheet_path):
            with open(sample_sheet_path, "r", encoding="utf-8") as f:
                lines = [line.strip() for line in f if line.strip()]
            in_section = False
            for line in lines:
                if line == "[DragenEnrichment_Data]":
                    in_section = True
                    continue
                if in_section:
                    if line == "Sample_ID":
                        continue
                    if line.startswith("[") and line.endswith("]"):
                        break
                    sample_ids.append(line)
        samplesize = len(sample_ids)
        return sequencing_date, sample_ids, samplesize
    return None, [], 0

def sync_created_chip(chipid):
    """
    當 chip 狀態為「晶片已創建」時，再次檢查是否已有 sequencing 資料；
    若有，更新 chip_info 並將 sample 匯入 sample_info（避免重複）。
    """
    sequencingDate, samplelist, samplesize = get_sequencing_info(chipid)

    # 沒有新資料就不動
    if not samplelist or not samplesize:
        return

    # 1) 更新 chip_info（僅限當前為「晶片已創建」）
    #    若你要更嚴謹，可另外加上判斷 sequencingDate 有值才更新
    sqlexe(
        """
        UPDATE chip_info
           SET "sequencingDate" = %s,
               "sampleSize"   = %s,
               "status"       = %s,
               "copycomplete" = %s
         WHERE "chipSNo"      = %s
           AND "status"       = '晶片已創建'
        """,
        [sequencingDate, sampleSize, "準備分析", True, chipid]
    )

    # 2) 匯入 sample_info（避免重複）
    # 需要資料庫上有唯一鍵 ("sampleid","chipSNo")
    for sample in samplelist:
        sqlexe(
            """
            INSERT INTO sample_info ("sampleid","chipSNo","status")
            VALUES (%s,%s,%s)
            ON CONFLICT ("sampleid","chipSNo") DO NOTHING
            """,
            [sample, chipid, "待分析"]
        )

def check_chip_status():
    sql = f"""select * from "chip_info" """
    chipdf = sqlquery(sql)

    # ---------- (A) 先處理「晶片已創建」→ 檢查是否已有 sequencing 資料可匯入 ----------
    created_df = chipdf[chipdf["status"] == "晶片已創建"]
    for _, row in created_df.iterrows():
        chipid = row["chipSNo"]
        try:
            sync_created_chip(chipid)
        except Exception as e:
            pass
    
    # ---------- (B) 再處理「檢測分析進行中」→ 根據 trace 決定後續 ----------
    in_progress_df = chipdf[chipdf['status'] == '檢測分析進行中']
    for index, row in in_progress_df.iterrows():
        chipid=row['chipSNo']
        analysisfolder=row['analysisfolder']
        tmp_analysis_path = os.path.join(settings.BASE_DIR, "wes", "static", "tmp", "analysis", analysisfolder)
        output_path = os.path.join(settings.BASE_DIR, "wes", "static", "analysis", analysisfolder)
        if not os.path.isdir(tmp_analysis_path) and not os.path.isdir(output_path):
            update_chip_info(chipid, '準備分析')
            continue
        if os.path.isdir(tmp_analysis_path):
            continue
        if not os.path.isdir(output_path):
            continue
        trace_path = os.path.join(output_path,"trace.txt")
        if trace_path and os.path.exists(trace_path):
            tracedf = pd.read_csv(trace_path, sep="\t")
            status_list = tracedf['status'].tolist()
            if len(status_list)!=0:
                all_completed = all(status == 'COMPLETED' for status in status_list)
            else:
                all_completed = False
        else:
            all_completed = False
        if all_completed: 
            #更新chipid狀態
            update_chip_info(chipid, '檢測分析已完成')
            process_result_summary(output_path, chipid)
            import_to_database(output_path)
        # else:
        #     #更新chipid狀態
        #     update_chip_info(chipid, '準備分析')
        #     # 發送郵件通知分析失敗
        #     # sending_mail_fail(chipid, output_path)

def pad(text):
    return text + (16 - len(text) % 16) * ' '

def encrypt(plain_text):
    """ AES 加密（ECB 模式） """
    cipher = AES.new(KEY, AES.MODE_ECB)
    padded_text = pad(plain_text)
    encrypted_bytes = cipher.encrypt(padded_text.encode())
    return base64.b64encode(encrypted_bytes).decode()

def decrypt(encrypted_text):
    """ AES 解密 """
    cipher = AES.new(KEY, AES.MODE_ECB)
    encrypted_bytes = base64.b64decode(encrypted_text)
    return cipher.decrypt(encrypted_bytes).decode().strip()

def safe_decrypt(x):
    try:
        if pd.notnull(x) and x != "":
            decoded_bytes = base64.b64decode(x)
            decrypted_value = decrypt(x)
            return decrypted_value
        else:
            return x
    except Exception as e:
        return x 


def safe_json_value(val):
    """
    處理存進josn中的NaN 
    2026/3/17 update
    """
    if isinstance(val, dict):
        return {k: safe_json_value(v) for k, v in val.items()}
    if isinstance(val, list):
        return [safe_json_value(v) for v in val]
    if pd.isna(val):
        return None
    if isinstance(val, (np.floating, float)):
        return float(val)
    if isinstance(val, (np.integer, int)):
        return int(val)
    return val

def variant_to_database(result_path,chipID):
    """
    Insert variant into psql
    2026/3/17 update
    """
    sql = f"""SELECT "sampleSNo" from sample_info where "chipSNo" =  '{chipID}' """
    samples=sqlquery(sql)
    for sample in samples['sampleSNo']:
        ## 讀取opencravat註解結果 (sqlite)
        print(f"load variant for sample {sample}")
        if(os.path.exists(f"{result_path}/{sample}")):
            conn = sqlite3.connect(f"{result_path}/{sample}/opencravat/{sample}.hard-filtered.sqlite")
            variant_table=pd.read_sql_query("SELECT * FROM variant ;", conn)
            ## 讀取inhouse-filteration的結果 (xlsx)
            filtered_table=pd.read_excel(f"{result_path}/{sample}/{sample}_filtered.xlsx")
            ## 合併兩表並以report欄位標示是否為篩選結果
            variant_table=variant_table.merge(filtered_table.assign(report=True),
                        left_on=['base__chrom','base__pos','base__ref_base','base__alt_base'],
                        right_on=['Chrom','Position','Ref Base','Alt Base'],
                        how='left')
            variant_table['report']= variant_table['report'].astype('boolean').fillna(False)
            variant_table['base__exonno']=variant_table['base__exonno'].astype("Int64") 
            ## insert rows to psql
            for i in range(0,variant_table.shape[0]):
                uniID=f"{sample}_{chipID}"
                variantID=i
                uniVID=f"{uniID}_{variantID}"
                chrom=variant_table['base__chrom'][i]
                pos=int(variant_table['base__pos'][i])
                ref=variant_table['base__ref_base'][i]
                alt=variant_table['base__alt_base'][i]
                dbsnp=variant_table['dbsnp__rsid'][i]
                report=bool(variant_table['report'][i])
                consequence=safe_json_value({
                    "gene":variant_table['base__hugo'][i],
                    "transcript":variant_table['base__transcript'][i],
                    "hgvsc":variant_table['base__cchange'][i],
                    "hgvsp":variant_table['base__achange'][i],
                    "exon": safe_json_value(variant_table['base__exonno'][i]),
                    "sequence_ontology":variant_table['base__so'][i]
                }) 
                population={
                    "gnomad":{
                        "global":safe_json_value({
                            "AC":variant_table['gnomad4__ac'][i],
                            "AN":variant_table['gnomad4__an'][i],
                            "AF":variant_table['gnomad4__af'][i],
                            "Homo":variant_table['gnomad4__nhomalt'][i],
                        }),
                        "AFR":safe_json_value({
                            "AC":variant_table['gnomad4__ac_afr'][i],
                            "AN":variant_table['gnomad4__an_afr'][i],
                            "AF":variant_table['gnomad4__af_afr'][i],
                            "Homo":variant_table['gnomad4__nhomalt_afr'][i],
                        }),
                        "AMR":safe_json_value({
                            "AC":variant_table['gnomad4__ac_amr'][i],
                            "AN":variant_table['gnomad4__an_amr'][i],
                            "AF":variant_table['gnomad4__af_amr'][i],
                            "Homo":variant_table['gnomad4__nhomalt_amr'][i],
                        }),
                        "EAS":safe_json_value({
                            "AC":variant_table['gnomad4__ac_eas'][i],
                            "AN":variant_table['gnomad4__an_eas'][i],
                            "AF":variant_table['gnomad4__af_eas'][i],
                            "Homo":variant_table['gnomad4__nhomalt_eas'][i],
                        }),
                        "SAS":safe_json_value({
                            "AC":variant_table['gnomad4__ac_sas'][i],
                            "AN":variant_table['gnomad4__an_sas'][i],
                            "AF":variant_table['gnomad4__af_sas'][i],
                            "Homo":variant_table['gnomad4__nhomalt_sas'][i],
                        }),
                        "FIN":safe_json_value({
                            "AC":variant_table['gnomad4__ac_fin'][i],
                            "AN":variant_table['gnomad4__an_fin'][i],
                            "AF":variant_table['gnomad4__af_fin'][i],
                            "Homo":variant_table['gnomad4__nhomalt_fin'][i],
                        }),
                        "NFE":safe_json_value({
                            "AC":variant_table['gnomad4__ac_nfe'][i],
                            "AN":variant_table['gnomad4__an_nfe'][i],
                            "AF":variant_table['gnomad4__af_nfe'][i],
                            "Homo":variant_table['gnomad4__nhomalt_nfe'][i],
                        })
                    }
                }
                clinvar=safe_json_value({
                    "clininical_significance":variant_table['clinvar__sig'][i],
                    "review_status":variant_table['clinvar__rev_stat'][i],
                    "clinvar_id":variant_table['clinvar__id'][i],
                    "significance_detail":variant_table['clinvar__sig_conf'][i],
                })
                vcf_info=safe_json_value({
                    "total_reads":variant_table['vcfinfo__tot_reads'][i],
                    "alt_reads":variant_table['vcfinfo__alt_reads'][i],
                    "allele_fraction":variant_table['vcfinfo__af'][i],
                    "quality":variant_table['vcfinfo__phred'][i],
                    "zygosity":variant_table['vcfinfo__zygosity'][i],
                    "filter":variant_table['vcfinfo__filter'][i]
                })
                prediction={
                    "spliceai":safe_json_value({
                        "score":variant_table.loc[i,variant_table.columns[variant_table.columns.str.contains('spliceai__ds')]].max(),
                        "class":"Pathogenic" if (variant_table.loc[i,variant_table.columns[variant_table.columns.str.contains('spliceai__ds')]]>0.5).any() else None
                    }),
                    "dbscsnv_ada":safe_json_value({
                        "score":variant_table['dbscsnv__ada_score'][i],
                        "class":"Pathogenic" if variant_table['dbscsnv__ada_score'][i]>0.7 else None
                    }),
                    "revel":safe_json_value({
                        "score":variant_table['revel__rankscore'][i],
                        "class":f"Pathogenic_{variant_table['revel__pp3_pathogenic'][i]}" if pd.notna(variant_table['revel__pp3_pathogenic'][i]) else None
                    }),
                    "vest4":safe_json_value({
                        "score":variant_table['vest__score'][i],
                        "class":f"Pathogenic_{variant_table['vest__pp3_pathogenic'][i]}" if pd.notna(variant_table['vest__pp3_pathogenic'][i]) else None
                    })
                }
                sqlexe(
                    """
                    INSERT INTO sample_small_variant ("UniVID", "UniID", "variantID", "chrom", "pos","ref_base","alt_base","dbsnp","report","consequence","population","clinvar","vcf_info","prediction")
                    VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                    ON CONFLICT ("UniVID") DO NOTHING
                    """,
                    [uniVID, uniID, variantID, chrom, pos, ref, alt, dbsnp, report, json.dumps(consequence), json.dumps(population), 
                    json.dumps(clinvar),json.dumps(vcf_info),json.dumps(prediction)]
                )
            print(f"Completed!")
        else:
            print(f"no analytic result is found for sample {sample}")
            
    
        


            




#########################################API####################################################
## WESCoreAnalysis Section
### API-T01
def coreanalysis(request):
    """
    2026/3/18 update
    """
    # Get filter metadata
    try:
        userid = request.POST['userid']
        finditem = request.POST['finditem']
        findspecimen = request.POST['findspecimen']
        findsample = request.POST['findsample']
        findchip = request.POST['findchip']
        findstarttime = request.POST['findstarttime']
        findendtime = request.POST['findendtime']
        finderrortype = request.POST['finderrortype']
        finddetail = request.POST['finddetail']

        if userid is None:
            return JsonResponse({"Code":500, "Msg":"Error found on server"})

        # Datatable niptcoreanalysis
        sql = """select * from wescoreanalysis"""
        twes = sqlquery(sql)

        # Filter by item
        '''
        if finditem != '全部':
            tnipt = tnipt[tnipt['testItem'] == finditem]
        '''
        # Filter by specimen
        '''
        if len(findspecimen) != 0:
            tnipt = tnipt[tnipt['specimenNumber'] == findspecimen]
        '''
        # Filter by sample (用前9碼去對應)
        if len(findsample) != 0:
            twes=twes[twes['sampleSNo'].str.contains(findsample)]
            '''
            if findsample.startswith("PC"):
                # PC → 前 6 碼
                tnipt = tnipt[tnipt['sampleSNo'].str[:6] == findsample[:6]]
            # elif findsample.startswith("NC"):
            #     # NC → 前 7 碼
            #     tnipt = tnipt[tnipt['sampleSNo'].str[:9] == findsample[:9]]
            else:
                # 其他 → 前 9 碼
                tnipt = tnipt[tnipt['sampleSNo'].str[:9] == findsample[:9]]
            '''

        # Filter by chip
        if len(findchip) != 0:
            twes = twes[twes['chipSNo'] == findchip]

        # Filter by time
        if len(findstarttime) != 0:
            starttime = pd.to_datetime(findstarttime)
            tnipt = tnipt[tnipt['sequencingDate'] >= starttime]

        if len(findendtime) != 0:
            endtime = pd.to_datetime(findendtime)
            tnipt = tnipt[tnipt['sequencingDate'] <= endtime]

        # Filter by detail
        '''
        if finddetail != '全部':
            tnipt = tnipt[tnipt['testDetail'] == finddetail]
        '''
        # Filter by detail
        '''
        if finderrortype != '無':
            tnipt =  tnipt[tnipt['redrawBlood'] == finderrortype]
        '''
        ## output
        #tnipt["totalRds"] = tnipt["totalRds"].apply(lambda x: f'{x:,}')
        #tnipt["uniMapRds"] = tnipt["uniMapRds"].apply(lambda x: f'{x:,}')
        twes = twes.fillna("-")
        twes = twes.replace('nan', '-').replace('', '-')
        #tnipt["NIPTTestResults"] = tnipt["NIPTTestResults"].replace('-', '低風險')
        twes["NIPTTestResults"] = twes["report_variant_count"]
        twes = convertTime(twes, 'sequencingDate')
        twes = convertTime(twes, 'analysisTime')
        core_table = twes.to_json(orient='records')
        core_table = json.loads(core_table)
        return JsonResponse({"Code":200, "Msg":{'user_id':userid, 'core_table':core_table}})
    except:
        return JsonResponse({"Code":500, "Msg": "No data found"})

### API-T02
def coreanalysis_download(request):
    list_testidx = request.POST['list_testidx']
    corelist = json.loads(list_testidx)

    sql = """SELECT * FROM niptcoreanalysis"""
    coredf = sqlquery(sql)
    coredf = coredf.loc[:, ~coredf.columns.duplicated()]
    coredf = coredf.drop(columns=["UniID"])
    tnipt = coredf[coredf["testidx"].isin(corelist)]

    #下載頁面模板中的coreanalysis_selected_cols，轉換成中文
    sql = """select "coreanalysis_selected_cols","coreanalysis_cols" from download_cols"""
    tnipt = download_and_convert(sql,tnipt, ['coreanalysis_selected_cols', 'coreanalysis_cols'])

    ## Make file
    now = datetime.now()
    output_name = f"niptINFO_{now.strftime('%Y%m%d_%H%M%S')}.xlsx"
    outpath = os.path.join(settings.BASE_DIR, 'nipt/static/tmp/niptdata', output_name)
    filepath = os.path.join('/static/tmp/niptdata', output_name)
    tnipt.to_excel(outpath, index=False)

    return JsonResponse({"Code":200, "Msg":{"core_durl":filepath}})

### API-T03
def coreanalysis_detail(request):
    userid = request.POST['user_id']
    sampleid = request.POST['sampleSNo']
    chipid = request.POST['chipSNo']

    if userid is None:
        return JsonResponse({"Code":500, "Msg":"Error found on server"})
    
    T03sql = f"""SELECT chrom, pos, ref_base, alt_base, 
        consequence->>gene as gene,
        consequence->>hgvsc as hgvsc,
        consequence->>hgvsp as hgvsp,
        consequence->>transcript as transcript,
        consequence->>transcript as transcript,
                        from sample_small_variant where "UniID" = '{sampleid}_{chipid}'"""
    T03all= sqlquery(T03sql)
    T03all = T03all.apply(lambda x: float(x) if isinstance(x, Decimal) else x)

    # Sample info section
    ## section 1
    niptSample1 = T03all[['sampleSNo','pregnantName', "redrawBlood"]]
    niptSample1 = {k: v[0] for k, v in niptSample1.to_dict(orient='list').items()}

    ## section2
    niptSample2 = T03all[['sampleSNo','pregnantName', "bloodCollectionTime", "identifyNo", "pregnantAge", "testItem", "testDetail", "specimenNumber", "analysisTime"]]
    niptSample2 = niptSample2.apply(lambda x: float(x) if isinstance(x, Decimal) else x)
    niptSample2['identifyNo'] = niptSample2['identifyNo'].apply(safe_decrypt)
    niptSample2 = convertTime(niptSample2, 'bloodCollectionTime')
    niptSample2 = convertTime(niptSample2, 'analysisTime')
    niptSample2 = {k: v[0] for k, v in niptSample2.to_dict(orient='list').items()}

    ### Check other sample with the same specimen
    specimenNumber = niptSample2['specimenNumber']
    othersql = f"""select "sampleSNo", "chipSNo" from "niptcoreanalysis" where "specimenNumber" = '{specimenNumber}' AND "sampleSNo" != '{sampleid}' """
    othersamples = sqlquery(othersql)
    if othersamples.empty:
        othersamplesdf = []
    else:
        othersamples['UID'] = othersamples['sampleSNo'] + '_' + othersamples['chipSNo']
        othersamplesdf = othersamples['UID']
        othersamplesdf = othersamplesdf.tolist()
    niptSample1['otherSampleSNoList'] = othersamplesdf

    # Quality info section
    ## section 1: QC summary
    niptQC1 = T03all[['sampleSNo','chipSNo','originalQualityControl','qualified']]
    niptQC1 = {k: v[0] for k, v in niptQC1.to_dict(orient='list').items()}
    allpass = T03all['qualified'].values[0]

    ## section2: QC detail (lack of FCPercent)
    niptQC2 = T03all.drop(['sampleSNo','chipSNo'],axis=1)
    niptQC2 = {k: v[0] for k, v in niptQC2.to_dict(orient='list').items()}

    # Analysis info section
    ## section 1: ID
    niptR1 = T03all[['sampleSNo','chipSNo','status','qualified']]
    niptR1 = {k: v[0] for k, v in niptR1.to_dict(orient='list').items()}

    ## section 2: result
    niptR2 = T03all[['sampleSNo','chipSNo','pregnantName','analysisProcess',"redrawBlood", "originalQualityControl"]]
    niptR2 = {k: v[0] for k, v in niptR2.to_dict(orient='list').items()}
    
    ## section 3: T13, T18, T21 更新為all chr table (same as CNV info section-niptR3)

    ## section 4: result (according to testItem...)
    ## Check item
    testItem = T03all['testItem'].values[0]
    testDetail = T03all['testDetail'].values[0]

    if testItem == 'NIPTPLUS' and testDetail == "FFAK":
        s44sql = f"""
        SELECT 
            "NIPTTestResults", 
            "CNVin38" as "CNVTestResults_in",
            "CNVout38" as "CNVTestResults_out",
            "microdeletion_all",
            CASE 
                WHEN (
                    ("T13_18_21" IS NOT NULL AND "T13_18_21" != 'nan') OR
                    ("SCA" IS NOT NULL AND "SCA" != 'nan') OR
                    ("microdeletion" IS NOT NULL AND "microdeletion" != 'nan') OR
                    ("RAAin38" IS NOT NULL AND "RAAin38" != 'nan') OR
                    ("RAAout38" IS NOT NULL AND "RAAout38" != 'nan')
                ) THEN 
                    CASE 
                        WHEN qualified = 'Yes' THEN '高風險'
                        ELSE '高風險, ' || REGEXP_REPLACE(qualified, '^No\\((.*)\\)$', '\\1') || ' fail'
                    END
                ELSE 
                    CASE 
                        WHEN qualified = 'Yes' THEN '無異常'
                        ELSE REGEXP_REPLACE(qualified, '^No\\((.*)\\)$', '\\1') || ' fail'
                    END
            END AS "NIPTReminderInformation"
        FROM "niptcoreanalysis"
        WHERE "sampleSNo" = '{sampleid}'
        AND "chipSNo" = '{chipid}';
        """
    elif testItem == 'NIPTPLUS' and testDetail == "FFAW":
        s44sql = f"""
        SELECT 
            "NIPTTestResults", 
            "CNVin38" as "CNVTestResults_in",
            "CNVout38" as "CNVTestResults_out",
            "microdeletion_all", 
            CASE 
                WHEN (
                ("T13_18_21" IS NOT NULL AND "T13_18_21"!= 'nan') 
                OR ("SCA" IS NOT NULL AND "SCA"!= 'nan') 
                OR ("microdeletion" IS NOT NULL AND "microdeletion"!= 'nan') 
                OR ("RAAin38" IS NOT NULL AND "RAAin38"!= 'nan') 
                OR ("RAAout38" IS NOT NULL AND "RAAout38"!= 'nan') 
                OR ("CNVin38" IS NOT NULL AND "CNVin38"!= 'nan') 
                OR ("CNVout38" IS NOT NULL AND "CNVout38"!= 'nan')
                ) THEN 
                    CASE 
                        WHEN qualified = 'Yes' THEN '高風險'
                        ELSE '高風險, ' || REGEXP_REPLACE(qualified, '^No\\((.*)\\)$', '\\1') || ' fail'
                    END
                ELSE 
                    CASE 
                        WHEN qualified = 'Yes' THEN '無異常'
                        ELSE REGEXP_REPLACE(qualified, '^No\\((.*)\\)$', '\\1') || ' fail'
                    END
            END AS "NIPTReminderInformation"
        FROM "niptcoreanalysis"
        WHERE "sampleSNo" = '{sampleid}'
        AND "chipSNo" = '{chipid}';
        """
    else:
        s44sql = f"""
        SELECT 
            "NIPTTestResults",
            "CNVin38" as "CNVTestResults_in",
            "CNVout38" as "CNVTestResults_out",
            "microdeletion_all", 
            CASE 
                WHEN (
                ("T13_18_21" IS NOT NULL AND "T13_18_21"!= 'nan') 
                OR ("SCA" IS NOT NULL AND "SCA"!= 'nan') 
                OR ("RAAin38" IS NOT NULL AND "RAAin38"!= 'nan') 
                OR ("RAAout38" IS NOT NULL AND "RAAout38"!= 'nan') 
                ) THEN 
                    CASE 
                        WHEN qualified = 'Yes' THEN '高風險'
                        ELSE '高風險, ' || REGEXP_REPLACE(qualified, '^No\\((.*)\\)$', '\\1') || ' fail'
                    END
                ELSE 
                    CASE 
                        WHEN qualified = 'Yes' THEN '無異常'
                        ELSE REGEXP_REPLACE(qualified, '^No\\((.*)\\)$', '\\1') || ' fail'
                    END
            END AS "NIPTReminderInformation"
        FROM "niptcoreanalysis"
        WHERE "sampleSNo" = '{sampleid}'
        AND "chipSNo" = '{chipid}';
        """

    niptR4 = sqlquery(s44sql)
    niptR4 = niptR4.replace('nan', '-').replace('', '-')
    niptR4 = niptR4.replace('-', '低風險')
    niptR4 = {k: v[0] for k, v in niptR4.to_dict(orient='list').items()}

    ## section 5: reference
    niptR5 = T03all[['RawReads','Total_rds','UniMap_rds','Unimap_gc','Duprate',"UniqPercent", "ff"]]
    niptR5 = {k: v[0] for k, v in niptR5.to_dict(orient='list').items()}

    # CNV info section
    ## Check CNV
    csql = f"""select * from "chip_info" where "chipSNo" = '{chipid}'"""
    chipdf = sqlquery(csql)
    analysisfolder= chipdf['analysisfolder'].iloc[0]

    sample_data = T03all
    result_path1 = os.path.join(settings.BASE_DIR, "nipt/static/tmp/analysis",analysisfolder,"result",sampleid,sampleid+".plots")
    result_path2 = os.path.join("static/tmp/analysis",analysisfolder,"result",sampleid,sampleid+".plots")
    chromosome_dict = {}
    niptR3 ={}
    #all Chr 
    microdeletion_list = sample_data['microdeletion_all'].apply(lambda x: x.split('\n')).tolist()
    for chr_label in list(range(1, 23)) + ['X', 'Y']:
        zscore = '-' if pd.isna(sample_data[f'{chr_label}_zscore'].values[0]) else round(float(sample_data[f'{chr_label}_zscore'].values[0]), 2)
        ratio = '-' if pd.isna(sample_data[f'{chr_label}_ratio'].values[0]) else round(float(sample_data[f'{chr_label}_ratio'].values[0])*100, 5)
        cnvplot = os.path.join(result_path1, f"chr{chr_label}.png")
        microdeletion_match = []
        for microdeletion in microdeletion_list:
            for entry in microdeletion:
                if f'loss({chr_label}:' in entry or f'gain({chr_label}:' in entry:
                    microdeletion_match.append(entry)
        if len(microdeletion_match) == 0:
            microdeletion_match.append('-')
        if os.path.exists(cnvplot):
            cnvplot=os.path.join(result_path2, f"chr{chr_label}.png")
            chromosome_dict[f'chr{chr_label}'] = {
                'zscore': zscore,
                'ratio': ratio,
                'CNVDownload': cnvplot,
                'microdeletion_all': '\n'.join(microdeletion_match)
            }
        niptR3[f'chr{chr_label}'] = {
            'zscore': zscore,
            'ratio': ratio
        }

    s7sql = f"""select * from cnvdata where "sampleSNo" = '{sampleid}' AND "chipSNo" = '{chipid}'"""
    cnvdata = sqlquery(s7sql)
    if cnvdata.shape[0] == 0:
        status = '-'
    else:
        status = ','.join(cnvdata['result'])

    zip_file = os.path.join(settings.BASE_DIR, "nipt/static/tmp/analysis",analysisfolder,"result",sampleid,sampleid+"_cnv.zip")
    zip_directory(result_path1, zip_file)
    zip_path=os.path.join("static/tmp/analysis",analysisfolder,"result",sampleid, f"{sampleid}_cnv.zip")

    niptR6 = {"sampleSNo":sampleid, "chipSNo":sample_data['chipSNo'].values[0],'qualified': allpass, "CNVTestResults":status, "CNVDownload":zip_path}
    final = {"sample_info":{"section1":niptSample1, "section2": niptSample2}, 
             "quality_info":{"section1":niptQC1, "section2":niptQC2}, 
             "analysis_info":{"section1":niptR1, "section2":niptR2, "section3":niptR3, "section4":niptR4, "section5":niptR5}, 
             "cnv_info":{"section1":niptR6,"section2":chromosome_dict}}        

    return JsonResponse({"Code":200, "Msg":final})

## Management Section
### API-C01
def chipsearch(request):
    try:
        check_chip_status()
        userid = request.POST['userid']
        findchip = request.POST['findchipSNo']
        findstarttime = request.POST['findstarttime']
        findendtime = request.POST['findendtime']
        findstatus = request.POST['findstatus']

        if userid is None:
            return JsonResponse({"Code":500, "Msg":"Error found on server"})

        # Datatable chips
        sql = """select * from chip_info"""
        chipdf = sqlquery(sql)
        # filter
        has_chip = len(findchip) != 0
        has_start_end = len(findstarttime) != 0 and len(findendtime) != 0
        has_status = len(findstatus) != 0 and findstatus != "全部"
        if has_chip:
            chipdf = chipdf[chipdf['chipSNo'] == findchip]
        if has_status:
            chipdf = chipdf[chipdf['status'] == findstatus]
        now = datetime.now()
        if has_start_end:
            starttime = pd.to_datetime(findstarttime)
            endtime = pd.to_datetime(findendtime)
            chipdf = chipdf[(chipdf['sequencingDate'] >= starttime) & (chipdf['sequencingDate'] <= endtime)]
        elif not has_chip and not has_status:
            # no filter using halfyear
            default_starttime = now - timedelta(days=180)
            default_endtime = now
            chipdf = chipdf[(chipdf['sequencingDate'].isnull()) |
                            ((chipdf['sequencingDate'] >= default_starttime) &
                             (chipdf['sequencingDate'] <= default_endtime))]

        ## output
        chipdf = chipdf.fillna("")
        chipdf = convertTime(chipdf, 'sequencingDate')
        core_table = chipdf.to_json(orient='records')
        core_table = json.loads(core_table)
        return JsonResponse({"Code":200, "Msg":{'user_id':userid, 'core_table': core_table}})
    
    except:
        return JsonResponse({"Code":500, "Msg": "No data found"})

### API-C02
def chip_add(request):
    chipid = request.POST['chipSNo_new'].strip()
    check_sql = 'SELECT 1 FROM chip_info WHERE "chipSNo" = %s LIMIT 1'
    check_result = sqlquery(check_sql, [chipid])
    if check_result is not None and not check_result.empty:
        return JsonResponse({"Code": 500, "Msg": f"已有相同晶片 {chipid} 於資料庫內"})
    sequencingDate, samplelist, samplesize = get_sequencing_info(chipid)

    try:
        # 沒有 sample 資料
        if not samplelist or samplesize == 0:
            sqlexe(
                """
                INSERT INTO chip_info ("chipSNo","sequencingDate","sampleSize","status")
                VALUES (%s,%s,%s,%s)
                """,
                [chipid, sequencingDate, samplesize, "晶片已創建"]
            )
            return JsonResponse({"Code": 200,"Msg": f"晶片 {chipid} 已建立，但沒有找到 sequencing 資料"})

        # 有 sample 資料
        sqlexe(
            """
            INSERT INTO chip_info ("chipSNo","sequencingDate","sampleSize","status","copycomplete")
            VALUES (%s,%s,%s,%s,%s)
            """,
            [chipid, sequencingDate, samplesize, "準備分析", True]
        )

        for sample in samplelist:
            UniID = f"{sample}_{chipid}"
            sqlexe(
                """
                INSERT INTO sample_info ("UniID","sampleSNo","chipSNo","status")
                VALUES (%s,%s,%s,%s)
                """,
                [UniID, sample, chipid, "待分析"]
            )

        return JsonResponse({"Code": 200, "Msg": f"新的晶片 {chipid} 已建立，共 {samplesize} 個 samples"})
    except Exception as e:
        return JsonResponse({"Code": 500, "Msg": "新增晶片失敗", "Error": str(e)})

### API-C03
def chip_edit(request):
    ori_chipid = request.POST['chipSNo']
    chipid = request.POST['chipSNo_new']
    try:
        sqlexe(f"""UPDATE chip_info SET "chipSNo" = %s where "chipSNo"= %s""", [chipid, ori_chipid])
        sequencingDate = get_sequencing_date_from_folders(chipid)
        sqlexe(f"""UPDATE chip_info SET "sequencingDate" = %s where "chipSNo"= %s""", [sequencingDate, chipid])
        return JsonResponse({"Code":200, "Msg":f"""{ori_chipid}晶片已更新名稱成{chipid}，頁面請重新整理!"""})
    except:
        return JsonResponse({"Code":500, "Msg":"晶片編輯失敗"})

### API-C04
def chip_delete(request):
    chip_list = json.loads(request.POST.get('list_chipSNo', '[]'))
    try:
        sqlexe(f"""DELETE FROM chip_info WHERE "chipSNo" IN %s""", [tuple(chip_list)])
        sqlexe(f"""DELETE FROM sample_info WHERE "chipSNo" IN %s""", [tuple(chip_list)])
        sqlexe(f"""DELETE FROM "Mapping" where "chipSNo"= %s""", [tuple(chip_list)])
        sqlexe(f"""DELETE FROM "QCanalysis" where "chipSNo"= %s""", [tuple(chip_list)])
        sqlexe(f"""DELETE FROM "niptdata" where "chipSNo"= %s""", [tuple(chip_list)])
        sqlexe(f"""DELETE FROM "cnvdata" where "chipSNo"= %s""", [tuple(chip_list)])
        sqlexe(f"""DELETE FROM "microdeletion" where "chipSNo"= %s""", [tuple(chip_list)])
        chip_str = ", ".join(chip_list)
        return JsonResponse({"Code":200, "Msg":f"""{chip_str}晶片與相關樣本紀錄皆從資歷料庫已刪除，頁面請重新整理！"""})
    except:
        return JsonResponse({"Code":500, "Msg":"並無此紀錄於資料庫內"})

### API-C05
def chip_renew(request):
    chipid = request.POST['chipSNo']
    try:
        sqlexe(f"""UPDATE chip_info SET "status" = '晶片已創建' where "chipSNo"= %s""", [chipid])
        #刪除所有結果
        sqlexe(f"""DELETE FROM sample_info where "chipSNo"= %s""", [chipid])
        sqlexe(f"""DELETE FROM "Mapping" where "chipSNo"= %s""", [chipid])
        sqlexe(f"""DELETE FROM "QCanalysis" where "chipSNo"= %s""", [chipid])
        sqlexe(f"""DELETE FROM "niptdata" where "chipSNo"= %s""", [chipid])
        sqlexe(f"""DELETE FROM "cnvdata" where "chipSNo"= %s""", [chipid])
        sqlexe(f"""DELETE FROM "microdeletion" where "chipSNo"= %s""", [chipid])
        # 確認晶片的樣本數量
        csql = f"""select * from "sample_info" where "chipSNo" = '{chipid}'"""
        chipdf = sqlquery(csql)
        sample_N = len(chipdf)
        # Update chipid status
        chipsql = f"""UPDATE "chip_info" SET "sampleSize" = %s WHERE "chipSNo" = %s"""
        sqlexe(chipsql, [sample_N, chipid])

        return JsonResponse({"Code":200, "Msg":{"status":"晶片已創建", "msg":f"""{chipid}晶片編號已回歸至創建狀態"""}})
    except:
        return JsonResponse({"Code":500, "Msg":"不明原因出錯，請洽管理人員!"})

### API-C06
def chip_cancel(request):
    chipid = request.POST['chipSNo']
    try:
        sqlexe(f"""UPDATE chip_info SET "status" = '晶片已作廢' where "chipSNo"= %s""", [chipid])
        # #更新sample_metadata supplement狀態-->樣品錄入
        # update_sample_metadata(chipid, '樣本錄入')
        sqlexe(f"""DELETE FROM sample_info where "chipSNo"= %s""", [chipid])
        sqlexe(f"""DELETE FROM "Mapping" where "chipSNo"= %s""", [chipid])
        sqlexe(f"""DELETE FROM "QCanalysis" where "chipSNo"= %s""", [chipid])
        sqlexe(f"""DELETE FROM "niptdata" where "chipSNo"= %s""", [chipid])
        sqlexe(f"""DELETE FROM "cnvdata" where "chipSNo"= %s""", [chipid])
        sqlexe(f"""DELETE FROM "microdeletion" where "chipSNo"= %s""", [chipid])
        # 確認晶片的樣本數量
        csql = f"""select * from "sample_info" where "chipSNo" = '{chipid}'"""
        chipdf = sqlquery(csql)
        sample_N = len(chipdf)
        # Update chipid status
        chipsql = f"""UPDATE "chip_info" SET "sampleSize" = %s WHERE "chipSNo" = %s"""
        sqlexe(chipsql, [sample_N, chipid])
        return JsonResponse({"Code":200, "Msg":{"status":"晶片已作廢","msg":f"""{chipid}晶片編號已作廢"""}})
    except:
        return JsonResponse({"Code":500, "Msg":"不明原因發生，請洽管理人員!"})

### API-C07
def chip_download(request):
    userid = request.POST['userid']
    chiplist_str = request.POST['list_chipSNo']
    chiplist = json.loads(chiplist_str)

    if userid is None:
        return JsonResponse({"Code":500, "Msg":"Error found on server"})

    # Datatable chips
    sql = """select "chipSNo", "sequencingDate", "sampleSize" , "status", "analysisTime" from chip_info"""
    chipdf = sqlquery(sql)
    sql = """select * from run_qc"""
    qcdf = sqlquery(sql)
    alldf=pd.merge(chipdf,qcdf,on="chipSNo",how='left')
    twes = alldf[alldf["chipSNo"].isin(chiplist)]

    params = (tuple(chiplist),)
    sql = """
    SELECT
        s."chipSNo",
        s."sampleSNo",
        s."status",
        q."sex",
        q."input_reads",
        q."10x",
        q."mean_target_coverage",
        q."uniformity_of_coverage",
        q."aligned_reads",
        q."aligned",
        q."enrichment",
        q."padded_enrichment"
    FROM "sample_info" AS s
    JOIN "QCanalysis" AS q
      ON q."UniID" = s."UniID"
    WHERE s."chipSNo" IN %s
    ORDER BY s."chipSNo" ASC, s."testidx" ASC
    """
    sampleqc = sqlquery(sql, params)


    ## Make file
    now = datetime.now()
    output_name = f"chipINFO_{now.strftime('%Y%m%d_%H%M%S')}.xlsx"
    outpath = os.path.join(settings.BASE_DIR, 'wes/static/tmp/chip', output_name)
    filepath = os.path.join('/static/tmp/chip', output_name)

    with pd.ExcelWriter(outpath, engine="openpyxl") as writer:
        twes.to_excel(writer, sheet_name="Chips", index=False)
        sampleqc.to_excel(writer, sheet_name="SampleQC", index=False)

    return JsonResponse({"Code":200, "Msg":{"core_durl":filepath}})

### API-C08
def chip_sampleDetail(request):
    chipid = request.POST['chipSNo']

    # Check chipid exist in database:
    csql = f"""select * from "chip_info" where "chipSNo" = '{chipid}'"""
    chipdf = sqlquery(csql)
    if chipdf.shape[0] != 1:
        return JsonResponse({"Code":500, "Msg":"並已無此晶片紀錄，請重新整理頁面或者通知管理人員!"})
    
    sql = """
    SELECT
        s."testidx",
        s."sampleSNo",
        s."chipSNo",
        s."status",
        q."sex",
        q."input_reads",
        q."10x",
        q."mean_target_coverage",
        q."uniformity_of_coverage",
        q."aligned_reads",
        q."aligned",
        q."enrichment",
        q."padded_enrichment"
    FROM "sample_info" AS s
    JOIN "QCanalysis" AS q
    ON q."UniID" = s."UniID"   -- 以 UniID 對齊，最安全不會錯位
    WHERE s."chipSNo" = %s
    ORDER BY s."testidx" ASC;
    """
    dfs = sqlquery(sql, [chipid])
    dfs = dfs.fillna("")
    dftmp = dfs.to_json(orient='records')
    dfout = json.loads(dftmp)

    return JsonResponse({"Code":200, "Msg":{"chipSNo": chipid, "core_table":dfout}})

### API-C09
def chip_qc(request):
    chipid = request.POST['chipSNo']

    # Check chipid exist in database:
    csql = f"""select * from "chip_info" where "chipSNo" = '{chipid}'"""
    chipdf = sqlquery(csql)
    if chipdf.shape[0] != 1:
        return JsonResponse({"Code":500, "Msg":"並已無此晶片QC紀錄，請做完分析後再查詢!"})
    analysisfolder= chipdf['analysisfolder'].iloc[0]

    csql = f"""select * from "run_qc" where "chipSNo" = '{chipid}'"""
    chipdf = sqlquery(csql)
    row = chipdf.iloc[0].to_dict()
    row.pop("chipSNo", None)
    row.pop("created_at", None)
    original_qc = {k: coerce(v) for k, v in row.items()}

    core_table = {
        "originalQC": original_qc,
    }
    return JsonResponse({"Code":200, "Msg":{"chipSNo": chipid,"core_table":core_table}})

## Analysis Management Section
### API-A01
def analysis_search(request):
    # print(os.listdir(NS2000_path))
    try:
        userid = request.POST['userid']
        findchip = request.POST['findchipSNo']
        findstarttime = request.POST['findcomputerstarttime']
        findendtime = request.POST['findcomputerendtime']
        findstatus = request.POST['findstatus']

        if userid is None:
            return JsonResponse({"Code":500, "Msg":"Error found on server"})

        # Datatable chips
        # sql = f"""select "testidx", "chipSNo", "sequencingDate", "sampleSize", "status" from "chip_info" """
        sql = f"""select * from "chip_info" """
        chipdf = sqlquery(sql)
        #check nextflow finished
        check_chip_status()

        # Filter by chip
        if len(findchip) != 0:
            chipdf = chipdf[chipdf['chipSNo'] == findchip]

        # Filter by time
        if len(findstarttime) !=0:
            starttime = pd.to_datetime(findstarttime)
            chipdf = chipdf[chipdf['sequencingDate'] >= starttime]

        if len(findendtime) != 0:
            endtime = pd.to_datetime(findendtime)
            chipdf = chipdf[chipdf['sequencingDate'] <= endtime]

        if len(findstatus) != 0:
            if findstatus != "全部":
                chipdf = chipdf[chipdf['status'] == findstatus]

        for idx, row in chipdf.iterrows():
            if pd.isna(row['sequencingDate']):
                chipid = row['chipSNo']
                sequencingDate = get_sequencing_date_from_folders(chipid)
                if sequencingDate:
                    chipdf.at[idx, 'sequencingDate'] = sequencingDate
                    sqlexe(f"""UPDATE chip_info SET "sequencingDate" = %s where "chipSNo"= %s""", [sequencingDate, chipid])

        ## output
        chipdf = chipdf.fillna("")
        chipdf = convertTime(chipdf, 'sequencingDate')
        core_table = chipdf.to_json(orient='records')
        core_table = json.loads(core_table)

        return JsonResponse({"Code":200, "Msg":{'user_id':userid, 'core_table':core_table}})
    except:
        return JsonResponse({"Code":500, "Msg": "No data found"})
    
### API-A02
def analysis_sample(request):
    chipid = request.POST['chipSNo']

    # Check chipid exist in database:
    csql = f"""select * from "chip_info" where "chipSNo" = '{chipid}'"""
    chipdf = sqlquery(csql)
    if chipdf.shape[0] != 1:
        return JsonResponse({"Code":500, "Msg":"並已無此晶片紀錄，請重新整理頁面或者通知管理人員!"})

    sql = f"""select "testidx", "sampleSNo", "status" from "sample_info" where "chipSNo" = '{chipid}' """
    dfs = sqlquery(sql)
    dfs['operation'] = "刪除"
    dfs = dfs.fillna("")
    dfs = dfs.replace('nan', '-').replace('', '-')
    dfout = dfs.to_json(orient='records')
    dfsout = json.loads(dfout)

    return JsonResponse({"Code":200, "Msg":{"chipSNo": chipid, "core_table":dfsout}})

### API-A03
def analysis_submit(request):
    chipid = request.POST['chipSNo']
    # Check chipid exist in database:
    csql = f"""select * from "chip_info" where "chipSNo" = '{chipid}'"""
    chipdf = sqlquery(csql)
    if chipdf.shape[0] != 1:
        return JsonResponse({"Code":500, "Msg":"並已無此晶片紀錄，請重新整理頁面或者通知管理人員!"})
    todate = datetime.now().strftime("%Y%m%d_%H%M%S")
    # Check sample data is uploaded to the database
    if chipdf['status'].iloc[0] == "準備分析":
        # Check NS2000_path chipid位置: 
        directories = sorted([d for d in os.listdir(NS2000_path) if os.path.isdir(os.path.join(NS2000_path, d))], reverse=True)
        folders = [directory for directory in directories if chipid in directory]
        if len(folders)==1:
            folder=folders[0]
            batch_path = os.path.join(NS2000_path, folder)  

            #創建資料夾
            analysisfolder = f"{chipid}_{todate}"
            tmp_analysis_path = os.path.join(settings.BASE_DIR, "wes", "static", "tmp", "analysis")
            job_path = os.path.join(tmp_analysis_path, analysisfolder)
            os.makedirs(job_path)

            #從資料庫中抓取資料做sample_list.csv
            sql = f"""SELECT "sampleSNo" FROM "sample_info" WHERE "chipSNo" = '{chipid}' AND "status" = '待分析'"""
            dfs = sqlquery(sql)
            dfs.to_csv(os.path.join(job_path, "sample_list.csv"), index=False)

            #執行nextflow程式碼
            try:
                os.chdir(tmp_analysis_path)
                configfile = os.path.join(settings.BASE_DIR,"bin/nextflow.config")
                date = datetime.now().strftime("%Y/%m/%d")
                nf = os.path.join(settings.BASE_DIR, "bin/wes_pipeline.nf")
                command = f"{nextflow} run -ansi-log false {nf} --input {batch_path} --output {job_path} -c {configfile} > {job_path}/std.out 2> {job_path}/std.err -bg -with-trace {job_path}/trace.txt"
                run_nextflow = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                time.sleep(10)
                
                # #chip_info 更新
                status="檢測分析進行中"
                sqlexe(f"""UPDATE chip_info SET "status" = %s, "analysisTime" = %s, "analysisfolder" = %s where "chipSNo"= %s""", [status, date, analysisfolder, chipid])

                return JsonResponse({"Code":200, "Msg": {"msg": f"{chipid}晶片已開始分析!!!!"}})
            except:
                return JsonResponse({"Code":500, "Msg":"分析執行有誤，請確認"})

        elif len(folders)==0:
            return JsonResponse({"Code":500, "Msg":f"NS2000 path 未有{chipid}的下機資料"})
    else:
        return JsonResponse({"Code":500, "Msg":"無法執行該晶片分析"})

### API-A07        
def analysis_deletesample(request):
    chipid = request.POST['chipSNo']
    sampleid = request.POST['sampleSNo']

    # Check chipid exist in database:
    csql = f"""select * from "chip_info" where "chipSNo" = '{chipid}'"""
    chipdf = sqlquery(csql)
    if chipdf.shape[0] != 1:
        return JsonResponse({"Code":500, "Msg":"並已無此晶片紀錄，請重新整理頁面或者通知管理人員!"})

    # delete only uniid, not the sample id due to different chips
    uniid = sampleid + '_' + chipid
    try:
        sql = f"""UPDATE "sample_info" SET "status" = %s WHERE "UniID" = %s"""
        sqlexe(sql, ["已刪除", uniid])
        return JsonResponse({"Code": 200, "Msg": f"{chipid}晶片內{sampleid}狀態為已刪除!"})
    except:
        return JsonResponse({"Code":500, "Msg":"並無此紀錄於資料庫內"})

## Management Section
### API-M01
def login(request):
    user_account = request.POST.get('user_account')
    user_pwd = request.POST.get('password')
    
    if user_account and user_pwd:
        user_account = user_account.lower().strip()
        # Check user_account whether is involved in wesusers
        sql = """SELECT "user_id", "user_account", "user_pwd", "role" FROM "wesusers" WHERE "user_account" = %s"""
        users = sqlquery(sql, [user_account])

        if not users.empty:
            readPWD = users['user_pwd'][0]
            check_password_result = check_password(user_pwd, readPWD)
            if check_password_result:
                user_id = int(users['user_id'][0])
                user_role = users['role'][0]
                now = datetime.now()
                try:
                    # Write login information to the log history
                    sqlexe("""INSERT INTO "weslog" ("time", "user_account", "content") VALUES (%s, %s, %s)""",
                           [now, user_account, "check-in"])
                    return JsonResponse({"Code": 200, "Msg": {"user_account": user_account, "user_id": user_id, "role": user_role}})
                except Exception as e:
                    return JsonResponse({"Code": 500, "Msg": "登入記錄失敗，請稍後重試！", "Error": str(e)})
            else:
                return JsonResponse({"Code": 500, "Msg": "登入失敗，請確認密碼是否正確!"})
        else:
            return JsonResponse({"Code": 500, "Msg": "登入失敗，請確認帳號是否正確!"})
    else:
        return JsonResponse({"Code": 500, "Msg": "帳號或密碼至少有一個沒有填入!"})

### API-M02
def userpage(request):
    user_id = request.POST['user_id']
    if not user_id:
        return JsonResponse({"Code":500, "Msg":"沒有帳號身分，請洽管理系統管理員"})
    
    sql = f"""SELECT "user_id", "user_account", "user_pwd", "role" from "wesusers" """
    users = sqlquery(sql)
    users.columns = ["id", "user_account", "user_pwd", "role"]
    users = users.fillna("")
    core_table = users.to_json(orient='records')
    core_table = json.loads(core_table)

    sql = f"""SELECT * from "wesrole" """
    roles = sqlquery(sql)
    roles = roles.fillna("")
    role_table = roles.to_json(orient='records')
    role_table = json.loads(role_table)

    return JsonResponse({"Code":200, "Msg":{"user_table":core_table,"role_list":role_table}})

### API-M03
def adduser(request):
    user_account = request.POST.get('user_account')
    user_pwd = request.POST.get('user_pwd')
    user_role = request.POST.get('role')
    user_account = user_account.lower().strip()
    # Check user_account exist in database
    sql = """SELECT "user_account", "user_id" FROM "wesusers" WHERE "user_account" = %s"""
    userdf = sqlquery(sql, [user_account])
    if userdf.shape[0] >= 1:
        return JsonResponse({"Code": 500, "Msg": "已有此帳號紀錄，無法新增"})  
    # Generate password hash
    hashed_pwd = make_password(user_pwd)
    try:
        sqlexe("""INSERT INTO "wesusers" ("user_account", "user_pwd", "role") VALUES (%s, %s, %s)""",
               [user_account, hashed_pwd, user_role])
        return JsonResponse({"Code": 200, "Msg": f"新的用戶{user_account}已新增於資料庫，頁面請重新整理!"})
    except Exception as e:
        return JsonResponse({"Code": 500, "Msg": "新增發生錯誤，請洽管理員!", "Error": str(e)})
    
### API-M04
def edituser(request):
    user_account = request.POST['user_account']
    user_id = request.POST['user_id']
    ischange_pwd = request.POST['ischange_pwd']
    user_pwd = request.POST['user_pwd']
    new_pwd = request.POST['new_pwd']
    user_role = request.POST['role']

    # Check if any data needs to be updated
    if ischange_pwd == '1' and new_pwd:
        # Generate new password hash
        hashed_new_pwd = make_password(new_pwd)
    else:
        hashed_new_pwd = None

    # Build update query
    update_query = """UPDATE "wesusers" SET """
    update_fields = []

    if hashed_new_pwd:
        update_fields.append(f'"user_pwd" = %s')
    if user_role is not None:
        update_fields.append(f'"role" = %s')

    if not update_fields:
        return JsonResponse({"Code": 500, "Msg": "無須更新資料!"})

    update_query += ', '.join(update_fields)
    update_query += f""" WHERE "user_account" = '{user_account}'"""

    # Prepare parameters for execution
    params = []
    if hashed_new_pwd:
        params.append(hashed_new_pwd)
    if user_role is not None:
        params.append(user_role)
    try:
        sqlexe(update_query, params)
        return JsonResponse({"Code": 200, "Msg": f"""帳號{user_account}中的資料已經更新於資料庫"""})
    except Exception as e:
        return JsonResponse({"Code": 500, "Msg": "系統發生錯誤，請洽管理人員!", "Error": str(e)})    

### API-M05
def deleteuser(request):
    user_account = request.POST['user_account']
    try:
        sqlexe(f"""DELETE FROM "wesusers" where "user_account"= %s""", [user_account])
        return JsonResponse({"Code":200, "Msg":f"""帳號{user_account}已從資料庫刪除！"""})
    except:
        return JsonResponse({"Code":500, "Msg":"並無此紀錄於資料庫內"})

### API-M06
def displayrole(request):
    user_id = request.POST['user_id']
    if not user_id:
        return JsonResponse({"Code":500, "Msg":"沒有帳號身分，請洽管理系統管理員"})
    
    sql = f"""SELECT * from "wesrole" """
    roles = sqlquery(sql)
    roles.columns = ['id','role','description']
    roles = roles.fillna("")
    role_table = roles.to_json(orient='records')
    role_table = json.loads(role_table)
    return JsonResponse({"Code":200, "Msg":role_table})

### API-M07
def displayversion(request):
    try:
        sql = f"""SELECT * from "wesversion" """
        version = sqlquery(sql)
        version = version.fillna("")
        version = convertTime(version, 'updateTime')
        version_table = version.to_json(orient='records')
        version_table = json.loads(version_table)
        return JsonResponse({"Code":200, "Msg":version_table})
    except:
        return JsonResponse({"Code":500, "Msg":"系統錯誤"})

### API-M08
def loginhistory(request):
    try:
        sql = f"""SELECT * from "weslog" WHERE time >= current_date - INTERVAL '1 month' """
        logindf = sqlquery(sql)
        logindf = convertTime(logindf, 'time')
        logindf = logindf.fillna("")
        login_table = logindf.to_json(orient='records')
        login_table = json.loads(login_table)
        return JsonResponse({"Code":200, "Msg":login_table})
    except:
        return JsonResponse({"Code":500, "Msg":"系統錯誤"})


### API-M09
def loginversion(request):
    try:
        sql = f"""SELECT * from "wesversion" ORDER BY "updateTime" DESC LIMIT 1"""
        versiondf = sqlquery(sql)
        vid = versiondf['vid'][0]
        return JsonResponse({"Code":200, "Msg":{"version":vid}})
    except:
        return JsonResponse({"Code":500, "Msg":"系統錯誤!"})

### API-M10
def systemlogout(request):
    user_account = request.POST['user_account']
    now=datetime.now()
    try:
        sqlexe(f"""INSERT INTO "weslog" ("time","user_account", "content") VALUES (%s, %s, %s)""", [now, user_account, "check-out"])
        return JsonResponse({"Code":200, "Msg":f"""你的帳號{user_account}已經從WES系統上登出了! 可以用新的帳號登入了!"""})
    except:
        return JsonResponse({"Code":500, "Msg":"登出失敗，不明原因"})

### API-M11
def addversion(request):
    version = request.POST['verSNo']
    content = request.POST['content']
    try:
        sqlexe(f"""INSERT INTO "wesversion" ("vid", "content") VALUES (%s, %s)""", [version, content])
        return JsonResponse({"Code":200, "Msg":f"""你的版號第{version}號已經更新到系統上了，頁面請重新整理!"""})
    except:
        return JsonResponse({"Code":500, "Msg":"更新失敗!"})
    
