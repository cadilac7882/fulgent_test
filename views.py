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
from django.db import transaction
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
import numpy as np
from typing import Optional, Tuple, Dict, Any, List
import re
from collections import defaultdict

#########################################setting####################################################
NS2000_path = '/bioinfo/NS2000/RD/' #NS2000_path chipid位置
annotators_path = "/home/bioinfo/ggawes/bin/modules/annotators"
nextflow="/opt/nextflow"
# fail_mail_list=['ChenGyiLin@BionetTX.com','KennethYang@GGA.ASIA']
# success_mail_list=['ChenGyiLin@BionetTX.com','KennethYang@GGA.ASIA']
fail_mail_list=['ChenGyiLin@BionetTX.com']
success_mail_list=['ChenGyiLin@BionetTX.com']
KEY = 'GGANIPTSYSTEM0123456789876543210'.encode()
#########################################function###################################################
def sqlexe(query, params=None, returnValue=False):
    """
    2026/03/20 update: accept arg"returnValue"
    """
    try:
        with connection.cursor() as cursor:
            cursor.execute(query, params)
            if returnValue:
                return_value=cursor.fetchone()[0]
            connection.commit()
        if returnValue:
            return(return_value)
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

def update_status(table: str, chipid: str, status: str):
    """
    Update status by chipSNo for chip_info or sample_info
    """
    sql = f'''UPDATE "{table}" SET "status" = %s WHERE "chipSNo" = %s'''
    sqlexe(sql, [status, chipid])

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
        # sample_sheet_path = os.path.join(folder_path, "SampleSheet.csv")
        matches = glob.glob(os.path.join(folder_path, "**", "SampleSheet.csv"), recursive=True)
        sample_sheet_path = matches[0] if matches else None
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
        [sequencingDate, samplesize, "準備分析", True, chipid]
    )

    # 2) 匯入 sample_info（避免重複）
    for sample in samplelist:
        UniID = f"{sample}_{chipid}"
        sqlexe(
            """
            INSERT INTO sample_info ("UniID","sampleSNo","chipSNo","status")
            VALUES (%s,%s,%s,%s)
            ON CONFLICT ("UniID") DO NOTHING
            """,
            [UniID, sample, chipid, "待分析"]
        )

def check_chip_status():
    '''
    3/25 update check trace.txt & database.txt 
    '''
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
    
    # ---------- (B) 再處理「檢測分析進行中」→ 根據 trace & database 決定後續 ----------
    in_progress_df = chipdf[chipdf['status'] == '檢測分析進行中']
    for index, row in in_progress_df.iterrows():
        chipid=row['chipSNo']
        analysisfolder=row['analysisfolder']
        tmp_analysis_path = os.path.join(settings.BASE_DIR, "wes", "static", "tmp", "analysis", analysisfolder)
        output_path = os.path.join(settings.BASE_DIR, "wes", "static", "analysis", analysisfolder)
        if not os.path.isdir(tmp_analysis_path) and not os.path.isdir(output_path):
            update_status("chip_info", chipid, '準備分析')
            continue
        if os.path.isdir(tmp_analysis_path):
            continue
        if not os.path.isdir(output_path):
            continue

        trace_path = os.path.join(output_path, "trace.txt")
        database_path = os.path.join(output_path, "database.txt")
        all_completed = False
        if os.path.exists(trace_path):
            tracedf = pd.read_csv(trace_path, sep="\t")
            status_list = tracedf["status"].tolist()
            if status_list:
                all_completed = all(status == "COMPLETED" for status in status_list)
        all_completed = all_completed and os.path.exists(database_path)

        # ---------- 更新狀態 ----------
        if all_completed:
            # 更新 chip / sample 狀態
            update_status("chip_info", chipid, "檢測分析已完成")
            update_status("sample_info", chipid, "分析完成")
        else:
            update_status("chip_info", chipid, "準備分析")
            # 發送郵件通知分析失敗
            sending_mail_fail(chipid, output_path)

    in_progress_df = chipdf[chipdf['status'] == '重跑分析進行中']
    for index, row in in_progress_df.iterrows():
        chipid=row['chipSNo']
        analysisfolder=row['analysisfolder']
        output_path = os.path.join(settings.BASE_DIR, "wes", "static", "analysis", analysisfolder)
        trace_path = os.path.join(output_path, "trace2.txt")
        all_completed = False
        has_fail = False
        if os.path.exists(trace_path):
            tracedf = pd.read_csv(trace_path, sep="\t")
            status_list = tracedf["status"].tolist()
            if status_list:
                completed_count = status_list.count("COMPLETED")
                all_completed = all(status == "COMPLETED" for status in status_list)
                has_fail = any(status == "FAILED" for status in status_list)
        # ---------- 更新狀態 ----------
        if all_completed:
            if completed_count==2:
                update_status("chip_info", chipid, "檢測分析已完成")
        elif has_fail:
            update_status("chip_info", chipid, "準備分析")
            sending_mail_fail(chipid, output_path)
        else:
            pass

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
    2026/3/20 update: deal with bool
    2026/3/25 update:  add nan & Decimal
    """
    if isinstance(val, dict):
        return {k: safe_json_value(v) for k, v in val.items()}
    if isinstance(val, list):
        return [safe_json_value(v) for v in val]
    if pd.isna(val):
        return None
    if isinstance(val, bool):
        return bool(val)
    if isinstance(val, (np.floating, float)):
        return float(val)
    if isinstance(val, (np.integer, int)):
        return int(val)
    if str(val).lower() == 'nan':
        return "-"
    if isinstance(val, Decimal):
        return float(val)
    return val

def batch_upsert(
        sql_query: str, 
        rows: list,
        batch_size: int = 1000):
    """
    Insert multiple rows batchly
    """
    if not rows:
        return
    with transaction.atomic():
            with connection.cursor() as cursor:
                for start in range(0, len(rows), batch_size):
                    end = start + batch_size
                    cursor.executemany(
                        sql_query, 
                        rows[start:end]
                    )

def variant_to_database(result_path,chipID):
    """
    Insert variant into psql
    2026/3/17 update
    2026/3/20 update: separate sample_small_variant to varaint, sample_variant and varirant_consequence,insert rows by batch_upsert function
    2026/3/23 update: alter insertion query for sample_variant and variant
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
            
            ## insert rows to database
            variant_cache = {}
            sample_variant_rows = []
            consequence_rows= []
            for _, row in variant_table.iterrows():
                row_dict=row.to_dict()
                row_dict=safe_json_value(row_dict)
                chrom=row_dict['base__chrom']
                pos=int(row_dict['base__pos'])
                ref_base=row_dict['base__ref_base']
                alt_base=row_dict['base__alt_base']

                key = (chrom, pos, ref_base, alt_base)
                

                if key not in variant_cache:
                    ## upsert rows to variant table and get variant_id
                    variant_id=sqlexe("""
                            INSERT INTO variant (chrom, pos, ref_base, alt_base)
                            VALUES (%s, %s, %s, %s)
                            ON CONFLICT (chrom, pos, ref_base, alt_base)
                            DO UPDATE SET chrom = EXCLUDED.chrom
                            RETURNING variant_id;
                            """, 
                            [chrom, pos, ref_base, alt_base],True)
                    variant_cache[key]=variant_id
                else:
                    variant_id=variant_cache[key]
                
                ## collect rows for sample_variant and consequence
                sample_id = f"{sample}_{chipID}"
                
                sample_variant_rows.append((
                   sample_id,variant_id,row_dict['vcfinfo__zygosity'],row_dict['vcfinfo__tot_reads'],row_dict['vcfinfo__alt_reads'],row_dict['vcfinfo__filter'],row_dict['vcfinfo__phred'],
                     row_dict['extra_vcf_info__FS'],row_dict['extra_vcf_info__QD'],row_dict['extra_vcf_info__SOR'],row_dict['extra_vcf_info__MQ'],row_dict['extra_vcf_info__MQRankSum'],row_dict['extra_vcf_info__ReadPosRankSum'],
                     bool(row_dict['report'])
                ))

                consequence_rows.append((
                    variant_id,row_dict['base__hugo'],row_dict['base__transcript'],row_dict['base__exonno'],row_dict['base__cchange'],row_dict['base__achange'],row_dict['base__so'],"ENSEMBL"
                ))

            ## batch upsert rows to sample_variant table    
            batch_upsert(
                """
                INSERT INTO sample_variant (
                    sample_id, variant_id, genotype, total_reads, alt_reads,
                    filter, quality, FS, QD, SOR, MQ, MQRankSum, ReadPosRankSum, report
                )
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                ON CONFLICT (sample_id, variant_id) DO UPDATE 
                SET
                    genotype       = EXCLUDED.genotype,
                    total_reads    = EXCLUDED.total_reads,
                    alt_reads      = EXCLUDED.alt_reads,
                    filter         = EXCLUDED.filter, 
                    quality        = EXCLUDED.quality, 
                    FS             = EXCLUDED.FS, 
                    QD             = EXCLUDED.QD, 
                    SOR            = EXCLUDED.SOR, 
                    MQ             = EXCLUDED.MQ, 
                    MQRankSum      = EXCLUDED.MQRankSum, 
                    ReadPosRankSum = EXCLUDED.ReadPosRankSum, 
                    report         = EXCLUDED.report;
                """,
                sample_variant_rows
            )

            batch_upsert(
                """
                INSERT INTO variant_consequence (
                    variant_id, gene, transcript, exon, hgvsc, hgvsp, impact, source
                )
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
                ON CONFLICT DO NOTHING;
                """,
                consequence_rows
            )
            print(f"Completed!")
        else:
            print(f"no analytic result is found for sample {sample}")

def cnv_to_database(result_path,chipID):
    """2026/3/30
    Insert cnv into psql
    """
    sql = f"""SELECT "sampleSNo" from sample_info where "chipSNo" =  '{chipID}' """
    samples=sqlquery(sql)
    for sample in samples['sampleSNo']:
        ## 讀取opencravat註解結果 (sqlite)
        print(f"load variant for sample {sample}")
        if(os.path.exists(f"{result_path}/{sample}")):
            variant_table=pd.read_table(f"{result_path}/{sample}/AnnotSV/{sample}.cnv.tsv",sep='\t')
            variant_table=variant_table[variant_table['Annotation_mode']=='split']
            variant_table['SV_chrom']=variant_table['SV_chrom'].astype(str)

            ## 讀取inhouse-filteration的結果 (xlsx)
            filtered_table=pd.read_excel(f"{result_path}/{sample}/{sample}_filtered.xlsx",sheet_name='CNV')
            filtered_table=filtered_table[['SV_chrom','SV_start','SV_end']]
            filtered_table['SV_chrom']=filtered_table['SV_chrom'].astype(str)

            ## 合併兩表並以report欄位標示是否為篩選結果
            variant_table=variant_table.merge(filtered_table.assign(report=True),
                        on=['SV_chrom','SV_start','SV_end'],
                        how='left')
            variant_table['report']= variant_table['report'].astype('boolean').fillna(False)
            
            ## 切割vcf information
            variant_table[['GT', 'SM', 'CN', 'BC', 'PE']] = (
                variant_table[variant_table['Samples_ID'].iloc[0]]
                .str.split(':', expand=True)
            )
            variant_table['CN'] = pd.to_numeric(variant_table['CN'], errors='coerce').astype('Int64')
            variant_table['BC'] = pd.to_numeric(variant_table['BC'], errors='coerce').astype('Int64')
            variant_table['SM'] = pd.to_numeric(variant_table['SM'], errors='coerce')
            variant_table['GT']=variant_table['GT'].map({
                '0/1':'het',
                '1/1':'hom',
                './1':'unknown'
            })

            ## combine trascript id
            variant_table['transcript']=variant_table['Tx'].astype(str) + '.' + variant_table['Tx_version'].astype(int).astype(str)

            ## 計算affect exons
            variant_table[['exons', 'affect_exons']] = variant_table.apply(
                lambda x: pd.Series(
                    count_exons(x['Location'], x['Exon_count']),
                    index=['exons', 'affect_exons']
                ),
                axis=1
            )

            ## 強制讓chrom中包含chr  
            variant_table['SV_chrom'] = variant_table['SV_chrom'].where(
                variant_table['SV_chrom'].str.startswith('chr'),
                'chr' + variant_table['SV_chrom']
            )

            ## insert rows to database
            variant_cache = {}
            sample_variant_rows = []
            consequence_rows= []
            for _, row in variant_table.iterrows():
                row_dict=row.to_dict()
                row_dict=safe_json_value(row_dict)
                chrom    =row_dict['SV_chrom']
                start_pos=int(row_dict['SV_start'])
                end_pos  =int(row_dict['SV_end'])
                cn_type  =row_dict['SV_type']

                key = (chrom, start_pos, end_pos, cn_type)
                

                if key not in variant_cache:
                    ## upsert rows to variant table and get variant_id
                    variant_id=sqlexe("""
                            INSERT INTO cnv (chrom, start_pos, end_pos, cn_type)
                            VALUES (%s, %s, %s, %s)
                            ON CONFLICT (chrom, start_pos, end_pos, cn_type)
                            DO UPDATE SET chrom = EXCLUDED.chrom
                            RETURNING variant_id;
                            """, 
                            [chrom, start_pos, end_pos, cn_type],True)
                    variant_cache[key]=variant_id
                else:
                    variant_id=variant_cache[key]
                
                ## collect rows for sample_variant and consequence
                sample_id = f"{sample}_{chipID}"
                
                sample_variant_rows.append((
                   sample_id,variant_id,row_dict['GT'],row_dict['FILTER'],row_dict['QUAL'],row_dict['CN'],row_dict['SM'],
                     row_dict['BC'],bool(row_dict['report']
                )))

                consequence_rows.append((
                    variant_id,row_dict['Gene_name'],row_dict['transcript'],row_dict['exons'],row_dict['affect_exons'],"RefSeq"
                ))

            ## batch upsert rows to sample_variant table    
            batch_upsert(
                """
                INSERT INTO sample_cnv (
                    sample_id, variant_id, genotype, filter, quality, cn, sm, bc, report
                )
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
                ON CONFLICT (sample_id, variant_id) DO UPDATE 
                SET
                    genotype       = EXCLUDED.genotype,
                    filter         = EXCLUDED.filter, 
                    quality        = EXCLUDED.quality, 
                    cn             = EXCLUDED.cn, 
                    sm             = EXCLUDED.sm, 
                    bc             = EXCLUDED.bc,
                    report         = EXCLUDED.report;
                """,
                sample_variant_rows
            )

            batch_upsert(
                """
                INSERT INTO cnv_consequence (
                    variant_id, gene, transcript, exons, affect_exons, source
                )
                VALUES (%s, %s, %s, %s, %s, %s)
                ON CONFLICT DO NOTHING;
                """,
                consequence_rows
            )
            print(f"Completed!")
        else:
            print(f"no analytic result is found for sample {sample}")

def count_exons(region, total_exons=None):
    """2026/3/30
    Count the number of exons overlapped by a genomic region.

    Parameters
    ----------
    region : str
        Region string, e.g. "exon3-exon5", "intron6-exon8", "txStart-exon3"
    total_exons : int, optional
        Total number of exons in the gene (required if txEnd is used)

    Returns
    -------
    tuple
        (number_of_exons, affected_exon_label)
    """

    region = region.lower()

    # txStart-txEnd => full gene
    if region == "txstart-txend":
        if total_exons is None:
            raise ValueError("total_exons must be provided for txStart-txEnd")
        return total_exons, "Full gene"

    left, right = region.split("-")

    def get_num(s):
        nums = re.sub(r"[^0-9]", "", s)
        return int(nums) if nums else None

    # left boundary
    if left.startswith("exon"):
        start_exon = get_num(left)
    elif left.startswith("intron"):
        start_exon = get_num(left) + 1
    elif left == "txstart":
        start_exon = 1
    else:
        raise ValueError(f"Invalid left boundary: {left}")

    # right boundary
    if right.startswith("exon"):
        end_exon = get_num(right)
    elif right.startswith("intron"):
        end_exon = get_num(right)
    elif right == "txend":
        if total_exons is None:
            raise ValueError("total_exons must be provided when using txEnd")
        end_exon = int(total_exons)
    else:
        raise ValueError(f"Invalid right boundary: {right}")

    # calculate
    if end_exon < start_exon:
        return 0, "No exon"
    elif end_exon == start_exon:
        return 1, f"Exon{start_exon}"
    elif total_exons is not None and (end_exon - start_exon + 1) == total_exons:
        return total_exons, "Full gene"
    else:
        return end_exon - start_exon + 1, f"Exon{start_exon}-{end_exon}"
            
####抓資料庫版本
TOP_TITLE_RE = re.compile(r'^title\s*:\s*(.+?)\s*$')
TOP_VERSION_RE = re.compile(r'^version\s*:\s*(.+?)\s*$')

def read_text(path: str) -> List[str]:
    """以 UTF-8 讀取檔案，回退 errors='replace'，回傳行清單。"""
    try:
        with open(path, "r", encoding="utf-8") as f:
            return f.read().splitlines()
    except UnicodeDecodeError:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            return f.read().splitlines()

def extract_top_level_title_version(path: str) -> Tuple[Optional[str], Optional[str]]:
    """
    只擷取「行首（無縮排）」的 title/version。
    - 跳過空行與以 # 開頭的註解行。
    - 去除包覆的單/雙引號。
    - 只取第一個匹配到的值。
    """
    title = None
    version = None
    for line in read_text(path):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        if line[:1].isspace():
            continue
        if title is None:
            m = TOP_TITLE_RE.match(line)
            if m:
                val = m.group(1).strip()
                if (val.startswith('"') and val.endswith('"')) or (val.startswith("'") and val.endswith("'")):
                    val = val[1:-1]
                title = val
        if version is None:
            m = TOP_VERSION_RE.match(line)
            if m:
                val = m.group(1).strip()
                if (val.startswith('"') and val.endswith('"')) or (val.startswith("'") and val.endswith("'")):
                    val = val[1:-1]
                version = val
        if title is not None and version is not None:
            break
    return title, version

def scan_yaml_files(base_dir: str):
    """
    base_dir 下所有 .yml/.yaml 檔，回傳結果清單。
    每筆包含：module(父資料夾名), file(檔名), title, version, path(完整路徑)
    """
    results = []
    for root, _, files in os.walk(base_dir):
        for fn in files:
            if not fn.lower().endswith((".yml", ".yaml")):
                continue
            full = os.path.join(root, fn)
            module = os.path.basename(os.path.dirname(full))
            title, version = extract_top_level_title_version(full)
            results.append({
                "title": title if title else "N/A",
                "version": version if version else "N/A",
            })
    results.sort(key=lambda r: (r["title"].lower()))
    return results

def get_latest_vid():
    sql = ''' SELECT vid FROM wesversion ORDER BY "updateTime" DESC LIMIT 1'''
    df = sqlquery(sql)
    if df is None or df.empty:
        return None
    return df.iloc[0]['vid']

def build_population(r):
    """
    輸入object array，編排 gnomad 的欄位
    """
    result = []
    GNOMAD_SUBPOPS = [
        ('African', 'afr'),
        ('American', 'amr'),
        ('East Asian', 'eas'),
        ('European Finnish', 'fin'),
        ('European Non Finnish', 'nfe'),
        ('South Asian', 'sas'),
        ('Total', ''),
    ]
    for name, suffix in GNOMAD_SUBPOPS:
        s = f"_{suffix}" if suffix else ""
        result.append({
            'subpopulation': name,
            'allele_frequency': r[f'gnomad4__af{s}'] if pd.notna(r[f'gnomad4__af{s}']) else '-',
            'allele_count': int(r[f'gnomad4__ac{s}']) if pd.notna(r[f'gnomad4__ac{s}']) else '-',
            'allele_number': int(r[f'gnomad4__an{s}']) if pd.notna(r[f'gnomad4__an{s}']) else '-',
            'allele_homo_alt': int(r[f'gnomad4__nhomalt{s}']) if pd.notna(r[f'gnomad4__nhomalt{s}']) else '-',
        })
    return result

def auto_detect_panel_name(sex,sample_id):
    '''
    deal with sex
    '''
    if sex:
        sample_sex = 'Male' if sex=='XY' else 'Female'
    else:
        sample_sex = 'Female'
    '''
    assign panel
    '''
    if sample_id.startswith('CSB'):
        return(f"CSB_{sample_sex}")
    elif sample_id.startswith('CSF'):
        return(f"CSF_{sample_sex}")
    elif sample_id.startswith('SA'):
        return(f"SA")
    else:
        return(f"CSF_Female")

def collapse_ps_events(df):
    df = df.copy()
    df["_ps_group"] = df["PS"].where(
        df["PS"].notna(),
        pd.Series(df.index, index=df.index)
    )
    rows = []
    for _, g in df.groupby("_ps_group"):
        row = {}
        row["gene"] = g["gene"].iloc[0]
        row["transcript"] = g["transcript"].iloc[0]
        row["genotype"] = g["genotype"].iloc[0]
        row["inheritance"] = g["inheritance"].iloc[0]
        row["condition_name"] = g["condition_name"].iloc[0]
        has_ps = g["PS"].notna().any()
        has_mutalyzer = g["mutalyzer_output"].notna().any()
        if has_ps and has_mutalyzer:
            row["variant_detail"] = (
                g["mutalyzer_output"].dropna().iloc[0]
                + ", "
                + g["protein_description"].dropna().iloc[0]
            )
        else:
            row["variant_detail"] = g["variant_detail"].iloc[0]
        rows.append(row)
    return pd.DataFrame(rows)

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
        findsample = request.POST['findsample']
        findchip = request.POST['findchipSNo']
        findstarttime = request.POST['findstarttime']
        findendtime = request.POST['findendtime']
        if userid is None:
            return JsonResponse({"Code":500, "Msg":"Error found on server"})

        # Datatable wescoreanalysis
        sql = """select "testidx","sampleSNo","chipSNo","sex","mean_target_coverage" AS "meanTargetCoverage", "10x" AS "pctOver10x", "QC_pass", "report_variant_count" AS "reportSmallVariantCount","report_cnv_count" AS "reportCNVCount", "sequencingDate", "analysisTime", "status" from wescoreanalysis"""
        twes = sqlquery(sql)

        # Filter by sample
        if len(findsample) != 0:
            twes=twes[twes['sampleSNo'].str.contains(findsample)]

        # Filter by chip
        if len(findchip) != 0:
            twes = twes[twes['chipSNo'] == findchip]

        # Filter by time
        if len(findstarttime) != 0:
            starttime = pd.to_datetime(findstarttime)
            twes = twes[twes['sequencingDate'] >= starttime]

        if len(findendtime) != 0:
            endtime = pd.to_datetime(findendtime)
            twes = twes[twes['sequencingDate'] <= endtime]

        ## output
        twes = twes.fillna("-")
        twes = twes.replace('nan', '-').replace('', '-')
        twes['originalQualityControl'] = twes['QC_pass'].apply(lambda x: "PASS" if x is True else "FAIL")
        twes = convertTime(twes, 'sequencingDate')
        twes = convertTime(twes, 'analysisTime')
        twes['reportSmallVariantCount']="SNV:"+twes['reportSmallVariantCount'].astype('str')+"\nCNV:"+twes['reportCNVCount'].astype('str')
        core_table = twes.to_json(orient='records')
        core_table = json.loads(core_table)
        #print(core_table)
        return JsonResponse({"Code":200, "Msg":{'user_id':userid, 'core_table':core_table}})
    except:
        return JsonResponse({"Code":500, "Msg": "No data found"})

### API-T02
def coreanalysis_download(request):
    list_testidx = request.POST['list_testidx']
    corelist = json.loads(list_testidx)
    sql = """SELECT * FROM wescoreanalysis"""
    coredf = sqlquery(sql)
    twes = coredf[coredf["testidx"].isin(corelist)]

    all_report_snv = []
    all_report_cnv = []
    for _, row in twes.iterrows():
        sampleSNo = row["sampleSNo"]
        chipSNo = row["chipSNo"]
        sample_id = f"{sampleSNo}_{chipSNo}"
        # ===== SNV =====
        sql_snv = """
            SELECT *
            FROM snv_annotation
            WHERE sample_id = %s
            AND report = TRUE
        """
        report_snv = sqlquery(sql_snv, [sample_id])
        if not report_snv.empty:
            report_snv["sample_id"] = sample_id
            all_report_snv.append(report_snv)
        # ===== CNV =====
        sql_cnv = """
            SELECT *
            FROM cnv_annotation
            WHERE sample_id = %s
            AND report = TRUE
            AND gene IN (SELECT gene FROM inheritance)
        """
        report_cnv = sqlquery(sql_cnv, [sample_id])
        if not report_cnv.empty:
            report_cnv["sample_id"] = sample_id
            all_report_cnv.append(report_cnv)

    if all_report_snv:
        report_snv_df = pd.concat(all_report_snv, ignore_index=True)
    else:
        report_snv_df = pd.DataFrame()

    if all_report_cnv:
        report_cnv_df = pd.concat(all_report_cnv, ignore_index=True)
    else:
        report_cnv_df = pd.DataFrame()

    ## Make file
    now = datetime.now()
    output_name = f"wesINFO_{now.strftime('%Y%m%d_%H%M%S')}.xlsx"
    outpath = os.path.join(settings.BASE_DIR, 'wes/static/tmp/wesdata', output_name)
    filepath = os.path.join('/static/tmp/wesdata', output_name)
    with pd.ExcelWriter(outpath, engine="openpyxl") as writer:
        report_snv_df.to_excel(writer,sheet_name="small_variant",index=False)
        report_cnv_df.to_excel(writer,sheet_name="CNV",index=False)

    return JsonResponse({"Code":200, "Msg":{"core_durl":filepath}})

### API-T03
def coreanalysis_detail(request):
    userid = request.POST['user_id']
    sampleid = request.POST['sampleSNo']
    chipid = request.POST['chipSNo']
    print([userid,sampleid,chipid])

    if userid is None:
        return JsonResponse({"Code":500, "Msg":"Error found on server"})

    sql = """
        SELECT
            "sampleSNo",
            "chipSNo",
            "sex",
            "mean_target_coverage",
            "10x",
            "QC_pass",
            "input_reads",
            "uniformity_of_coverage",
            "aligned_reads",
            "aligned",
            "enrichment",
            "padded_enrichment"
        FROM wescoreanalysis
        WHERE "sampleSNo" = %s
        AND "chipSNo"   = %s
    """
    twes = sqlquery(sql, [sampleid, chipid])
    twes['originalQualityControl'] = twes['QC_pass'].apply(lambda x: "PASS" if x is True else "FAIL")
         
    # Quality info section
    ## section 1: QC summary
    wesQC1 = twes[["sex" ,'mean_target_coverage','10x','originalQualityControl']]
    wesQC1 = {k: v[0] for k, v in wesQC1.to_dict(orient='list').items()}
    wesQC1 = safe_json_value(wesQC1)

    ## section2: QC detail
    wesQC2 = twes[['input_reads',"aligned_reads" ,'aligned','uniformity_of_coverage','enrichment','padded_enrichment']]
    wesQC2 = {k: v[0] for k, v in wesQC2.to_dict(orient='list').items()}
    wesQC2 = safe_json_value(wesQC2)

    # Analysis info section
    ## load opencravat annotation table
    analysisfolder=sqlquery(
        '''
        SELECT analysisfolder FROM chip_info WHERE "chipSNo" = %s 
        ''',
        [chipid])['analysisfolder'].loc[0]
    file_path = os.path.join(settings.BASE_DIR, "wes", "static", "analysis", analysisfolder)
    db_path = f"{file_path}/{sampleid}/opencravat/{sampleid}.hard-filtered.sqlite"

    try:
        with sqlite3.connect(db_path) as conn:
            snv_ann_table = pd.read_sql_query(
                "SELECT * FROM variant_update;",
                conn
            )
    except sqlite3.OperationalError as e:
        snv_ann_table = None
        print(f"[{sampleid}] SQLite error: {e}")
    except Exception as e:
        snv_ann_table = None
        print(f"[{sampleid}] Unexpected error: {e}")

    ## section 1: SNV reuslts
    sample_id = f"{sampleid}_{chipid}"
    sql = """
        SELECT *
        FROM snv_annotation
        WHERE sample_id = %s
            AND report = TRUE
     """
    report_snv = sqlquery(sql, [sample_id]) 
    # print(report_snv.columns)
    
    report_snv['hgvsp']=report_snv['hgvsp'].apply(lambda x: 'p.?' if pd.isna(x) else x)
    report_snv['variant_detail']=report_snv.apply(lambda x: f"{x['hgvsc']}, {x['hgvsp']}",axis=1)
    report_snv['genotype'] = report_snv['genotype'].replace({'unknown':'-','het':'Heterozygous','hom':'Homozygous'})
    # snv_report = report_snv[["gene" ,'mane_refseq_tx','variant_detail','genotype','inheritance','condition_name']]
    snv_report = report_snv.merge(snv_ann_table,
                       left_on=['chrom','pos','ref_base','alt_base'],
                       right_on=['base__chrom','base__pos','base__ref_base','base__alt_base'],
                       how='left')
    cols_keep = [
        "gene","mane_refseq_tx","variant_detail","genotype",
        "inheritance","condition_name",
        "PS","mutalyzer_output","protein_description"
    ]
    snv_report = snv_report.loc[:, cols_keep].copy()
    snv_report = snv_report.rename(columns={'mane_refseq_tx':'transcript'})
    snv_report = collapse_ps_events(snv_report) 
    snv_report = snv_report.to_dict(orient="records")

    ## section 2: CNV reuslts
    sql = """
        SELECT *
        FROM cnv_annotation
        WHERE sample_id = %s
            AND report = TRUE
            AND gene IN (SELECT gene FROM inheritance)
     """
    report_cnv = sqlquery(sql, [sample_id])
    report_cnv['variant_detail'] = report_cnv.apply(lambda x: f"{x['affect_exons']} {x['cn_type']}",axis=1)
    report_cnv['genotype'] = report_cnv['genotype'].replace({'unknown':'-','het':'Heterozygous','hom':'Homozygous'})
    cnv_report = report_cnv[["gene" ,'transcript','variant_detail','genotype','inheritance','condition_name']]
    cnv_report = cnv_report.to_dict(orient="records")

    # small variant section
    ## query variants inside gene panel
    sql='''
        SELECT *
        FROM snv_annotation
        WHERE sample_id = %s
            AND gene IN (SELECT gene FROM gene_panel WHERE panel_name = %s) 
        '''
    
    select_panel=auto_detect_panel_name(wesQC1['sex'],sample_id)
    
    snv_in_panel = sqlquery(sql, [sample_id,select_panel])
    snv_in_panel['location'] = snv_in_panel.apply(lambda x: f"{x['chrom']}:{str(x['pos'])}:{x['ref_base']}:{x['alt_base']}",axis=1)
    snv_in_panel['hgvsp'] = snv_in_panel['hgvsp'].apply(lambda x: 'p.?' if pd.isna(x) else x)
    snv_in_panel['variant_detail'] = snv_in_panel.apply(lambda x: f"{x['hgvsc']}, {x['hgvsp']}",axis=1)
    snv_in_panel['genotype'] = snv_in_panel['genotype'].replace({'unknown':'-','het':'Heterozygous','hom':'Homozygous'})

    ## load opencravat annotation table
    analysisfolder=sqlquery(
        '''
        SELECT analysisfolder FROM chip_info WHERE "chipSNo" = %s 
        ''',
        [chipid])['analysisfolder'].loc[0]
    file_path = os.path.join(settings.BASE_DIR, "wes", "static", "analysis", analysisfolder)
    db_path = f"{file_path}/{sampleid}/opencravat/{sampleid}.hard-filtered.sqlite"

    try:
        with sqlite3.connect(db_path) as conn:
            snv_ann_table = pd.read_sql_query(
                "SELECT * FROM variant;",
                conn
            )
    except sqlite3.OperationalError as e:
        snv_ann_table = None
        print(f"[{sampleid}] SQLite error: {e}")
    except Exception as e:
        snv_ann_table = None
        print(f"[{sampleid}] Unexpected error: {e}")

    ## merge two tables
    snv_in_panel = snv_in_panel.merge(snv_ann_table,
                       left_on=['chrom','pos','ref_base','alt_base'],
                       right_on=['base__chrom','base__pos','base__ref_base','base__alt_base'],
                       how='left')
    cols = ['clinvar__sig_conf','clinvar__sig','clinvar__rev_stat','clinvar__id']
    snv_in_panel[cols] = snv_in_panel[cols].fillna('-')
    snv_in_panel['allele_balance'] = snv_in_panel['vcfinfo__af'].round(3)

    ## assign info dict
    records = snv_in_panel.to_dict(orient='records')
    snv_in_panel['info'] = [
        {
            'coding_detail': {
                'total_reads': r['total_reads'],
                'allele_reads': r['alt_reads'],
                'allele_balance': r['allele_balance'],
                'genotype': r['genotype'],
                'filter': r['filter'],
                'quality': r['quality'],
                'fs': r['extra_vcf_info__FS'],
                'qd': r['extra_vcf_info__QD'],
                'sor': r['extra_vcf_info__SOR'],
                'mq': r['extra_vcf_info__MQ'],
                'mqranksum': r['extra_vcf_info__MQRankSum'],
                'readposranksum': r['extra_vcf_info__ReadPosRankSum'],
            },
            'consequence': {
                'gene': r['gene'],
                'transcript': r['mane_refseq_tx'],
                'exon': int(r['exon']) if pd.notna(r['exon']) else '-',
                'hgvsc': r['hgvsc'],
                'hgvsp': r['hgvsp'],
                'impact': r['impact'],
            },
            'clinvar': {
                'clinical_significance': r['clinvar__sig_conf'] if pd.notna(r['clinvar__sig_conf']) else r['clinvar__sig'],
                'review_status': r['clinvar__rev_stat'],
                'clinvar_id': r['clinvar__id'],
            },
            'population': build_population(r),
            'prediction':[
                {
                    'tool':'REVEL',
                    'score':r['revel__rankscore'] if pd.notna(r['revel__rankscore']) else '-',
                    'class':f"Pathogenic_{r['revel__pp3_pathogenic']}" if pd.notna(r['revel__pp3_pathogenic']) else '-'},
                {
                    'tool':'VEST4',
                    'score':r['vest__score'] if pd.notna(r['vest__score']) else '-',
                    'class':f"Pathogenic_{r['vest__pp3_pathogenic']}" if pd.notna(r['vest__pp3_pathogenic']) else '-'},
                {
                    'tool':'dbscsnv_ada',
                    'score':r['dbscsnv__ada_score'] if pd.notna(r['dbscsnv__ada_score']) else '-',
                    'class':"Pathogenic" if r['dbscsnv__ada_score']>0.7 else '-'},
                {
                    'tool':'spliceai',
                    'score':max(r['spliceai__ds_ag'],r['spliceai__ds_al'],r['spliceai__ds_dg'],r['spliceai__ds_dl']) if pd.notna(r['spliceai__ds_ag']) else '-',
                    'class':"Pathogenic" if max(r['spliceai__ds_ag'],r['spliceai__ds_al'],r['spliceai__ds_dg'],r['spliceai__ds_dl'])>0.5 else '-'
                },
            ]
        }
        for r in records
    ]

    snv_in_panel['chr']=snv_in_panel['chrom']
    ## section 1: reported variants
    report_variants_table = snv_in_panel[snv_in_panel['report']]
    report_variants_table = report_variants_table[['chr','location','gene','mane_refseq_tx','variant_detail','total_reads','genotype','fulgent_report_times','info']]
    report_variants_table = report_variants_table.rename(columns={'mane_refseq_tx':'transcript'})
    report_variants_table = report_variants_table.to_dict(orient="records")

    ## section 2: other variants within gene panel
    other_variants_table = snv_in_panel[~snv_in_panel['report']]
    other_variants_table = other_variants_table[['chr','location','gene','mane_refseq_tx','variant_detail','total_reads','genotype','fulgent_report_times','info']]
    other_variants_table = other_variants_table.rename(columns={'mane_refseq_tx':'transcript'})
    other_variants_table = other_variants_table.to_dict(orient="records")

    other_variants_table_by_chr = defaultdict(list)
    for var in other_variants_table:
        other_variants_table_by_chr[var["chr"]].append(var)

    # copy number variant section
    ## section1: reported cnv
    report_cnv['location']=report_cnv.apply(lambda x: f"{x['chrom']}:{str(x['start_pos'])}:{x['end_pos']}:{x['cn_type']}",axis=1)
    report_cnv_table = report_cnv[['location','gene','transcript','variant_detail','quality','genotype']]
    report_cnv_table = {k: v for k, v in report_cnv_table.to_dict(orient='list').items()}
    report_cnv_table = safe_json_value(report_cnv_table)

    ## section2: other cnv
    sql='''
        SELECT *
        FROM cnv_annotation
        WHERE sample_id = %s
            AND report = FALSE
            AND gene in (select gene from inheritance)
        '''
    other_cnv = sqlquery(sql, [sample_id])
    other_cnv['location']=other_cnv.apply(lambda x: f"{x['chrom']}:{str(x['start_pos'])}:{x['end_pos']}:{x['cn_type']}",axis=1)
    other_cnv['variant_detail']=other_cnv.apply(lambda x: f"{x['affect_exons']} {x['cn_type']}",axis=1)
    other_cnv_table = other_cnv[['location','gene','transcript','variant_detail','quality','genotype']]
    other_cnv_table = {k: v for k, v in other_cnv_table.to_dict(orient='list').items()}
    other_cnv_table = safe_json_value(other_cnv_table)

    final = {
        "quality_info":{"section1":wesQC1, "section2":wesQC2},
        "analysis_info":{"section1":snv_report, "section2":cnv_report},
        "snv_info":{"section1":report_variants_table,"section2":other_variants_table_by_chr},
        #"cnv_info":{"section1":[],"section2":[]}
    }
    print(f"report variants:{len(report_variants_table)}\nother variants:{len(other_variants_table)}\nused panel:{select_panel}") 
    #print(final)
    return JsonResponse({"Code":200, "Msg": final})

### API-T04
def change_variant_report_status(request):
    status_dict = {'add':True,'remove':False}
    try:
        sampleid = request.POST['sampleSNo']
        chipid = request.POST['chipSNo']
        location = request.POST['location']
        variants = json.loads(location)
        status = request.POST['status']
        variant_type = request.POST['type']
        print(type(location))
        print(sampleid, chipid, location, status, variant_type)
        
        success_transactions = []
        failed_transactions = []
        for variant in variants:
            if variant_type == 'snv':
                chrom, pos, ref, alt = variant.split(":")
                
                try:
                    variant_id = sqlexe(
                        """
                        SELECT variant_id FROM variant WHERE chrom=%s and pos=%s and ref_base=%s and alt_base=%s
                        """,
                        [chrom, int(pos), ref, alt],True)
                    try:
                        sqlexe(
                            """
                            UPDATE sample_variant SET report=%s WHERE sample_id=%s AND variant_id=%s
                            """,
                        [status_dict[status],f"{sampleid}_{chipid}",variant_id])
                        success_transactions.append({variant:'success'})
                    except Exception as e:
                        failed_transactions.append({variant:"Error when changing status"})
                except Exception as e:
                    failed_transactions.append({variant:"No variant id is fetched in database"})
                
            else:
                chrom, start_pos, end_pos, cn_type = variant.split(":")
                try:
                    variant_id = sqlexe(
                        """
                        SELECT variant_id FROM cnv WHERE chrom=%s and start_pos=%s and end_pos=%s and cn_type=%s
                        """,
                        [chrom, int(start_pos), int(end_pos), cn_type],True)
                    try:
                        sqlexe(
                            """
                            UPDATE sample_cnv SET report=%s WHERE sample_id=%s AND variant_id=%s
                            """,
                        [status_dict[status],f"{sampleid}_{chipid}",variant_id])
                        success_transactions.append({variant:'success'})
                    except Exception as e:
                        failed_transactions.append({variant:"Error when changing status"})

                except Exception as e:
                    failed_transactions.append({variant:"No variant id is fetched in database"})
        print(success_transactions)
        return JsonResponse({"Code":200, "Msg": f"{status} variants\nSuccess:{len(success_transactions)}\nFail:{len(failed_transactions)}"})
    except:
        return JsonResponse({"Code":500, "Msg": "change variant report status failed!"})

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

    vid = get_latest_vid()
    if vid is None:
        return JsonResponse({"Code": 500, "Msg": "wesversion vid not get"})

    try:
        # 沒有 sample 資料
        if not samplelist or samplesize == 0:
            sqlexe(
                """
                INSERT INTO chip_info ("chipSNo","sequencingDate","sampleSize","status","vid")
                VALUES (%s,%s,%s,%s,%s)
                ON CONFLICT ("chipSNo") DO NOTHING
                """,
                [chipid, sequencingDate, samplesize, "晶片已創建", vid]
            )
            return JsonResponse({"Code": 200,"Msg": f"晶片 {chipid} 已建立，但沒有找到 sequencing 資料"})

        # 有 sample 資料
        sqlexe(
            """
            INSERT INTO chip_info ("chipSNo","sequencingDate","sampleSize","status","copycomplete","vid")
            VALUES (%s,%s,%s,%s,%s,%s)
            ON CONFLICT ("chipSNo") DO NOTHING
            """,
            [chipid, sequencingDate, samplesize, "準備分析", True, vid]
        )

        for sample in samplelist:
            UniID = f"{sample}_{chipid}"
            sqlexe(
                """
                INSERT INTO sample_info ("UniID","sampleSNo","chipSNo","status")
                VALUES (%s,%s,%s,%s)
                ON CONFLICT ("UniID") DO NOTHING
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
        sqlexe(f"""DELETE FROM "run_qc" where "chipSNo"= %s""", [tuple(chip_list)])
        sqlexe(f"""DELETE FROM "QCanalysis" where "chipSNo"= %s""", [tuple(chip_list)])
        chip_str = ", ".join(chip_list)
        return JsonResponse({"Code":200, "Msg":f"""{chip_str}晶片與相關樣本紀錄皆從資歷料庫已刪除，頁面請重新整理！"""})
    except:
        return JsonResponse({"Code":500, "Msg":"並無此紀錄於資料庫內"})

### API-C05
def chip_renew(request):
    chipid = request.POST['chipSNo']
    try:
        sqlexe(f"""UPDATE chip_info SET "status" = '晶片已創建' where "chipSNo"= %s""", [chipid])
        # check_chip_status()
        return JsonResponse({"Code":200, "Msg":{"status":"晶片已創建", "msg":f"""{chipid}晶片編號已回歸至創建狀態"""}})
    except:
        return JsonResponse({"Code":500, "Msg":"不明原因出錯，請洽管理人員!"})

### API-C06
def chip_cancel(request):
    chipid = request.POST['chipSNo']
    try:
        sqlexe(f"""UPDATE chip_info SET "status" = '晶片已作廢' where "chipSNo"= %s""", [chipid])
        sqlexe(f"""DELETE FROM sample_info where "chipSNo"= %s""", [chipid])
        sqlexe(f"""DELETE FROM "QCanalysis" where "chipSNo"= %s""", [chipid])
        sqlexe(f"""DELETE FROM "run_qc" where "chipSNo"= %s""", [chipid])
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
    status = chipdf.iloc[0]["status"]
    if status == '檢測分析已完成':
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
    else:
        sql = """
        SELECT
            "testidx",
            "sampleSNo",
            "chipSNo",
            "status"
        FROM "sample_info"
        WHERE "chipSNo" = %s
        ORDER BY "testidx" ASC;
        """
    
    dfs = sqlquery(sql, [chipid])
    dfs = dfs.replace('nan', '-').replace('', '-')
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
    vsql = """SELECT vid FROM wesversion ORDER BY "updateTime" DESC LIMIT 1"""
    versiondf = sqlquery(vsql)
    latest_vid = versiondf.iloc[0]["vid"]
    if chipdf.shape[0] != 1:
        return JsonResponse({"Code":500, "Msg":"並已無此晶片紀錄，請重新整理頁面或者通知管理人員!"})
    todate = datetime.now().strftime("%Y%m%d_%H%M%S")
    status = chipdf['status'].iloc[0]
    vid = chipdf['vid'].iloc[0]
    aid = f"{chipid}-{vid}"
    chipanalysisfolder = chipdf['analysisfolder'].iloc[0]
    # Check sample data is uploaded to the database
    if status == "準備分析":
        # Check NS2000_path chipid位置: 
        directories = sorted([d for d in os.listdir(NS2000_path) if os.path.isdir(os.path.join(NS2000_path, d))], reverse=True)
        folders = [directory for directory in directories if chipid in directory]
        if len(folders)==1:
            folder=folders[0]
            batch_path = os.path.join(NS2000_path, folder)  

            if pd.isna(chipanalysisfolder):
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
                    sqlexe('UPDATE chip_info SET "status"=%s,"analysisTime"=%s,"analysisfolder"=%s WHERE "chipSNo"=%s',
                        [status, date, analysisfolder, chipid])
                    sqlexe('INSERT INTO chip_analysis ("aid","chipSNo","analysisTime","analysisfolder","vid") VALUES (%s,%s,%s,%s,%s)',
                        [aid, chipid, date, analysisfolder, vid])
                    return JsonResponse({"Code":200, "Msg": {"msg": f"{chipid}晶片已開始分析!!!!"}})
                except:
                    return JsonResponse({"Code":500, "Msg":"分析執行有誤，請確認"})
            else:
                if latest_vid != vid:
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
                        sqlexe('UPDATE chip_info SET "status"=%s,"analysisTime"=%s,"analysisfolder"=%s WHERE "chipSNo"=%s',
                            [status, date, analysisfolder, chipid])
                        sqlexe('INSERT INTO chip_analysis ("aid","chipSNo","analysisTime","analysisfolder","vid") VALUES (%s,%s,%s,%s,%s)',
                            [aid, chipid, date, analysisfolder, vid])
                        return JsonResponse({"Code":200, "Msg": {"msg": f"{chipid}晶片已開始分析新版本!!!!"}})
                    except:
                        return JsonResponse({"Code":500, "Msg":"新版本分析執行有誤，請確認"})
                else:
                    tmp_analysis_path = os.path.join(settings.BASE_DIR, "wes", "static", "tmp", "analysis")
                    analysis_path = os.path.join(settings.BASE_DIR, "wes", "static", "analysis")
                    job_path = os.path.join(analysis_path, chipanalysisfolder)

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
                        command = f"{nextflow} run -ansi-log false {nf} --input {batch_path} --output {job_path} --only_filter -c {configfile} > {job_path}/std2.out 2> {job_path}/std2.err -bg -with-trace {job_path}/trace2.txt"
                        run_nextflow = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                        time.sleep(10)
                        
                        # #chip_info 更新
                        status="重跑分析進行中"
                        sqlexe('UPDATE chip_info SET "status"=%s,"analysisTime"=%s WHERE "chipSNo"=%s',
                            [status, date, chipid])
                        sqlexe('UPDATE chip_analysis SET "analysisTime"=%s WHERE "aid"=%s',
                            [date, aid])

                        return JsonResponse({"Code":200, "Msg": {"msg": f"{chipid}晶片已開始重跑分析!!!!"}})
                    except:
                        return JsonResponse({"Code":500, "Msg":"重跑分析執行有誤，請確認"})
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
        items = scan_yaml_files(annotators_path)
        db_version = "; ".join([f'{it["title"]}-{it["version"]}' for it in items])
        if db_version:
            db_version += ";"
        sqlexe('INSERT INTO "wesversion" ("vid", "content", "db_version") VALUES (%s, %s, %s)',[version, content, db_version])
        return JsonResponse({"Code":200, "Msg":f"""你的版號第{version}號已經更新到系統上了，頁面請重新整理!"""})
    except:
        return JsonResponse({"Code":500, "Msg":"更新失敗!"})
    
### API-M12
def gene_panel(request):
    try:
        sql = f"""SELECT 
                    g.panel_name,
                    g.gene,
                    i.inheritance,
                    i.condition_name
                FROM gene_panel g
                LEFT JOIN inheritance i
                ON g.gene = i.gene
                ORDER BY g.panel_name, g.gene;
                """
        gene_panel_df = sqlquery(sql)

        panel_dict = defaultdict(list)
        for _, row in gene_panel_df.iterrows():
            panel_dict[row["panel_name"]].append({
                "gene": row["gene"],
                "inheritance": row["inheritance"],
                "condition_name": row["condition_name"]
            })

        panel_with_count = {}
        for panel_name, gene_list in panel_dict.items():
            panel_with_count[panel_name] = {
                "count": len(gene_list),
                "genes": gene_list
            }

        return JsonResponse({"Code":200, "Msg": panel_with_count})
    except:
        return JsonResponse({"Code":500, "Msg":"gene_panel error!"})
    
