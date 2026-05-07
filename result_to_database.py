import os
import sys
import pandas as pd
import numpy as np
import psycopg2
from datetime import datetime
import sqlite3
from psycopg2.extras import execute_batch
import re
import urllib.request as req
import requests
import time
import json

class ImportToPostgres:
    def __init__(self, db_info: dict):
        self.connection = psycopg2.connect(
            host=db_info['host'],
            port=db_info.get('port', 5432),
            database=db_info['database'],
            user=db_info['user'],
            password=db_info['password']
        )
        self.connection.autocommit = False

    def close(self):
        if self.connection:
            self.connection.close()

    def main_flow(self, output_path: str):
        chipID = os.path.basename(output_path).split("_")[0]
        self.import_run_qc(output_path)
        sex_dict = self.sample_qc_to_database(output_path)
        self.variant_to_database(output_path, chipID)
        self.cnv_to_database(output_path, chipID, sex_dict)
        self.write_done_file(output_path)
        self.variant_acmg_check

    # ----------------------------
    # DB connection
    # ----------------------------
    def sqlexe(self, query, params=None, fetchone=False):
        try:
            with self.connection.cursor() as cursor:
                cursor.execute(query, params)
                result = cursor.fetchone() if fetchone else None
            self.connection.commit()
            return result
        except Exception as e:
            self.connection.rollback()
            raise RuntimeError(f"Database error: {e}")

    def sqlquery(self, query, params=None):
        try:
            with self.connection.cursor() as cursor:
                cursor.execute(query, params)
                columns = [col[0] for col in cursor.description]
                rows = cursor.fetchall()
                return pd.DataFrame(rows, columns=columns)
        except Exception as e:
            pass
            return pd.DataFrame()

    def batch_upsert(self, sql_query, rows, batch_size=1000):
        if not rows:
            return
        try:
            with self.connection.cursor() as cursor:
                for i in range(0, len(rows), batch_size):
                    cursor.executemany(
                        sql_query,
                        rows[i:i + batch_size]
                    )
            self.connection.commit()
        except Exception as e:
            self.connection.rollback()
            raise
       
    # ----------------------------
    # import to database
    # ----------------------------
    def import_run_qc(self, result_path: str):
        print(f"import run qc....")
        csv_path = os.path.join(result_path, "run_qc.csv")
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"run_qc.csv not found: {csv_path}")

        df = pd.read_csv(csv_path, dtype=str)
        rename_map = {
            'chipid': 'chipSNo',
            'Average_%>=Q30': 'average_q30',
            'Yield_Gbp': 'yield_gbp',
            '%Clusters_PF': 'clusters_ph',
            '%Occupied': 'occupied'
        }
        df.rename(columns=rename_map, inplace=True)
        row = df.iloc[0]
        sql = """
            INSERT INTO run_qc
                ("chipSNo", "average_q30", "clusters_ph", "occupied", "yield_gbp")
            VALUES (%s, %s, %s, %s, %s)
            ON CONFLICT ("chipSNo") DO NOTHING
        """

        self.sqlexe(
            sql,
            [
                row.get('chipSNo'),
                row.get('average_q30'),
                row.get('clusters_ph'),
                row.get('occupied'),
                row.get('yield_gbp')
            ]
        )

    def sample_qc_to_database(self, result_path):
        print(f"import sample qc....")
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
        sex_dict = {}
        for index, row in sampleqc_df.iterrows():
            data = tuple(map(str, row))
            self.sqlexe(sql, data)
            sex_dict[row["sampleSNo"]] = row["sex"]
        return sex_dict

    def variant_to_database_old(self, result_path, chipID):
        """
        Insert variant into psql
        2026/3/17 update
        2026/3/20 update: separate sample_small_variant to varaint, sample_variant and varirant_consequence,insert rows by batch_upsert function
        2026/3/23 update: alter insertion query for sample_variant and variant
        """
        sql = f"""SELECT "sampleSNo" from sample_info where "chipSNo" =  '{chipID}' """
        samples = self.sqlquery(sql)
        for sample in samples['sampleSNo']:
            ## 讀取opencravat註解結果 (sqlite)
            print(f"load variant for sample {sample}")
            if(os.path.exists(f"{result_path}/{sample}")):
                conn = sqlite3.connect(f"{result_path}/{sample}/opencravat/{sample}.hard-filtered.sqlite")
                variant_table=pd.read_sql_query("SELECT * FROM variant_update;", conn)
                # variant_table=pd.read_sql_query("SELECT * FROM variant ;", conn)
                # ## 讀取inhouse-filteration的結果 (xlsx)
                # filtered_table=pd.read_excel(f"{result_path}/{sample}/{sample}_filtered.xlsx")
                # ## 合併兩表並以report欄位標示是否為篩選結果
                # variant_table=variant_table.merge(filtered_table.assign(report=True),
                #             left_on=['base__chrom','base__pos','base__ref_base','base__alt_base'],
                #             right_on=['Chrom','Position','Ref Base','Alt Base'],
                #             how='left')
                # variant_table['report']= variant_table['report'].astype('boolean').fillna(False)
                variant_table['base__exonno']=variant_table['base__exonno'].astype("Int64") 
                
                ## insert rows to database
                variant_cache = {}
                sample_variant_rows = []
                consequence_rows= []
                for _, row in variant_table.iterrows():
                    row_dict=row.to_dict()
                    row_dict=self.safe_json_value(row_dict)
                    chrom=row_dict['base__chrom']
                    pos=int(row_dict['base__pos'])
                    ref_base=row_dict['base__ref_base']
                    alt_base=row_dict['base__alt_base']

                    key = (chrom, pos, ref_base, alt_base)

                    if key not in variant_cache:
                        ## upsert rows to variant table and get variant_id
                        variant_id = self.sqlexe("""
                                INSERT INTO variant (chrom, pos, ref_base, alt_base)
                                VALUES (%s, %s, %s, %s)
                                ON CONFLICT (chrom, pos, ref_base, alt_base)
                                DO UPDATE SET chrom = EXCLUDED.chrom
                                RETURNING variant_id;
                                """, 
                                [chrom, pos, ref_base, alt_base],True)
                        variant_cache[key]=variant_id
                    
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
                self.batch_upsert(
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

                self.batch_upsert(
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

    def variant_to_database(self, result_path, chipID):
        """
        Insert variant into psql
        LOGIC 100% IDENTICAL to variant_to_database_old
        - schema unchanged
        - SQL semantics unchanged
        Optimized only by:
        * limiting SQLite columns
        * limiting Excel columns
        * using execute_batch
        """

        sql = """SELECT "sampleSNo" FROM sample_info WHERE "chipSNo" = %s"""
        samples = self.sqlquery(sql, [chipID])
        print(chipID)
        if samples.empty:
            print("⚠️ No samples found for this chipID")
            return

        for sample in samples["sampleSNo"]:
            print(f"load variant for sample {sample}")

            sample_dir = os.path.join(result_path, sample)
            if not os.path.exists(sample_dir):
                print(f"no analytic result is found for sample {sample}")
                continue

            # ========== SQLite: ONLY required columns ==========
            sqlite_path = f"{sample_dir}/opencravat/{sample}.hard-filtered.sqlite"
            # conn = sqlite3.connect(sqlite_path)
            conn = sqlite3.connect(
                f"file:{sqlite_path}?mode=ro&immutable=1",
                uri=True
            )
            variant_table = pd.read_sql_query(
                """
                SELECT
                    base__chrom,
                    base__pos,
                    base__ref_base,
                    base__alt_base,
                    dbsnp__rsid,
                    vcfinfo__zygosity,
                    vcfinfo__tot_reads,
                    vcfinfo__alt_reads,
                    vcfinfo__filter,
                    vcfinfo__phred,
                    extra_vcf_info__FS,
                    extra_vcf_info__QD,
                    extra_vcf_info__SOR,
                    extra_vcf_info__MQ,
                    extra_vcf_info__MQRankSum,
                    extra_vcf_info__ReadPosRankSum,
                    base__hugo,
                    base__transcript,
                    base__exonno,
                    base__cchange,
                    base__achange,
                    base__so,
                    report
                FROM variant_update;
                """,
                conn
            )
            conn.close()

            # # ========== Excel: ONLY key + report ==========
            # filtered_table = pd.read_excel(
            #     f"{sample_dir}/{sample}_filtered.xlsx",
            #     usecols=["Chrom", "Position", "Ref Base", "Alt Base"]
            # )
            # filtered_table["report"] = True

            # # ========== merge ==========
            # variant_table = variant_table.merge(
            #     filtered_table,
            #     left_on=["base__chrom", "base__pos", "base__ref_base", "base__alt_base"],
            #     right_on=["Chrom", "Position", "Ref Base", "Alt Base"],
            #     how="left"
            # )

            variant_table["report"] = (
                variant_table["report"]
                .astype("boolean")
                .fillna(False)
            )
            variant_table["base__exonno"] = variant_table["base__exonno"].astype("Int64")

            records = variant_table.to_dict("records")

            variant_cache = {}
            sample_variant_rows = []
            consequence_rows = []

            sample_id = f"{sample}_{chipID}"

            # ========== DB transaction per sample ==========
            try:
                with self.connection.cursor() as cursor:
                    for r in records:
                        r = self.safe_json_value(r)

                        chrom = r["base__chrom"]
                        pos = int(r["base__pos"])
                        ref_base = r["base__ref_base"]
                        alt_base = r["base__alt_base"]
                        rsid = r["dbsnp__rsid"]

                        key = (chrom, pos, ref_base, alt_base)

                        # --- SAME LOGIC AS ORIGINAL ---
                        if key not in variant_cache:
                            cursor.execute(
                                """
                                INSERT INTO variant (chrom, pos, ref_base, alt_base, rsid)
                                VALUES (%s, %s, %s, %s, %s)
                                ON CONFLICT (chrom, pos, ref_base, alt_base)
                                DO UPDATE SET
                                    rsid = COALESCE(variant.rsid, EXCLUDED.rsid)
                                RETURNING variant_id;
                                """,
                                (chrom, pos, ref_base, alt_base, rsid)
                            )
                            variant_id = cursor.fetchone()[0]
                            # cursor.execute(
                            #     """
                            #     INSERT INTO variant (chrom, pos, ref_base, alt_base)
                            #     VALUES (%s, %s, %s, %s)
                            #     ON CONFLICT (chrom, pos, ref_base, alt_base)
                            #     DO UPDATE SET chrom = EXCLUDED.chrom
                            #     RETURNING variant_id;
                            #     """,
                            #     (chrom, pos, ref_base, alt_base)
                            # )
                            # variant_id = cursor.fetchone()[0]
                            variant_cache[key] = variant_id
                        else:
                            variant_id = variant_cache[key]

                        sample_variant_rows.append((
                            sample_id,
                            variant_id,
                            r["vcfinfo__zygosity"],
                            r["vcfinfo__tot_reads"],
                            r["vcfinfo__alt_reads"],
                            r["vcfinfo__filter"],
                            r["vcfinfo__phred"],
                            r["extra_vcf_info__FS"],
                            r["extra_vcf_info__QD"],
                            r["extra_vcf_info__SOR"],
                            r["extra_vcf_info__MQ"],
                            r["extra_vcf_info__MQRankSum"],
                            r["extra_vcf_info__ReadPosRankSum"],
                            bool(r["report"])
                        ))

                        consequence_rows.append((
                            variant_id,
                            r["base__hugo"],
                            r["base__transcript"],
                            r["base__exonno"],
                            r["base__cchange"],
                            r["base__achange"],
                            r["base__so"],
                            "ENSEMBL"
                        ))

                    # ========== batch inserts ==========
                    execute_batch(
                        cursor,
                        """
                        INSERT INTO sample_variant (
                            sample_id, variant_id, genotype, total_reads, alt_reads,
                            filter, quality, FS, QD, SOR, MQ, MQRankSum, ReadPosRankSum, report
                        )
                        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                        ON CONFLICT (sample_id, variant_id) DO UPDATE SET
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
                        sample_variant_rows,
                        page_size=1000
                    )

                    execute_batch(
                        cursor,
                        """
                        INSERT INTO variant_consequence (
                            variant_id, gene, transcript, exon, hgvsc, hgvsp, impact, source
                        )
                        VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
                        ON CONFLICT DO NOTHING;
                        """,
                        consequence_rows,
                        page_size=1000
                    )

                self.connection.commit()
                print(f"Completed variant import for sample {sample}")

            except Exception as e:
                self.connection.rollback()
                raise RuntimeError(f"Variant import failed for sample {sample}: {e}")
    
    def gt_to_genotype(self, gt, chrom, sex):
        if sex == 'XY' and chrom == 'X':
            if gt in ('1', '1/1'):
                return 'hemi'
            else:
                return 'unknown'
        if gt == '1/1':
            return 'hom'
        elif gt == '0/1':
            return 'het'
        return 'unknown'

    def cnv_to_database(self, result_path, chipID, sex_dict):
        """2026/3/30
        Insert cnv into psql
        """
        sql = f"""SELECT "sampleSNo" from sample_info where "chipSNo" =  '{chipID}' """
        samples = self.sqlquery(sql)
        for sample in samples['sampleSNo']:
            ## 讀取opencravat註解結果 (sqlite)
            print(f"load cnv for sample {sample}")
            if(os.path.exists(f"{result_path}/{sample}")):
                variant_table=pd.read_table(f"{result_path}/{sample}/AnnotSV/{sample}.cnv_update.tsv",sep='\t')
                variant_table=variant_table[variant_table['Annotation_mode']=='split']
                variant_table['SV_chrom']=variant_table['SV_chrom'].astype(str)

                # ## 讀取inhouse-filteration的結果 (xlsx)
                # filtered_table=pd.read_excel(f"{result_path}/{sample}/{sample}_filtered.xlsx",sheet_name='CNV')
                # filtered_table=filtered_table[['SV_chrom','SV_start','SV_end']]
                # filtered_table['SV_chrom']=filtered_table['SV_chrom'].astype(str)

                # ## 合併兩表並以report欄位標示是否為篩選結果
                # variant_table=variant_table.merge(filtered_table.assign(report=True),
                #             on=['SV_chrom','SV_start','SV_end'],
                #             how='left')
                variant_table['report']= variant_table['report'].astype('boolean').fillna(False)
                
                ## 切割vcf information
                variant_table[['GT', 'SM', 'CN', 'BC', 'PE']] = (
                    variant_table[variant_table['Samples_ID'].iloc[0]]
                    .str.split(':', expand=True)
                )
                variant_table['CN'] = pd.to_numeric(variant_table['CN'], errors='coerce').astype('Int64')
                variant_table['BC'] = pd.to_numeric(variant_table['BC'], errors='coerce').astype('Int64')
                variant_table['SM'] = pd.to_numeric(variant_table['SM'], errors='coerce')
                sex = sex_dict.get(sample)
                variant_table['GT'] = variant_table.apply(
                    lambda r: self.gt_to_genotype(
                        gt=r['GT'],
                        chrom=r['SV_chrom'],
                        sex=sex
                    ),
                    axis=1
                )
                # variant_table['GT'] = (
                #     variant_table['GT']
                #     .map({
                #         '0/1': 'het',
                #         '1/1': 'hom',
                #         './1': 'unknown',
                #         '1': 'hom'
                #     })
                #     .fillna('unknown')
                # )

                ## combine trascript id
                variant_table['transcript']=variant_table['Tx'].astype(str) + '.' + variant_table['Tx_version'].astype(int).astype(str)

                ## 計算affect exons
                variant_table[['exons', 'affect_exons']] = variant_table.apply(
                    lambda x: pd.Series(
                        self.count_exons(x['Location'], x['Exon_count']),
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
                    row_dict = row.to_dict()
                    row_dict = self.safe_json_value(row_dict)
                    chrom    = row_dict['SV_chrom']
                    start_pos= int(row_dict['SV_start'])
                    end_pos  = int(row_dict['SV_end'])
                    cn_type  = row_dict['SV_type']

                    key = (chrom, start_pos, end_pos, cn_type)
                    

                    if key not in variant_cache:
                        ## upsert rows to variant table and get variant_id
                        variant_id = self.sqlexe("""
                                INSERT INTO cnv (chrom, start_pos, end_pos, cn_type)
                                VALUES (%s, %s, %s, %s)
                                ON CONFLICT (chrom, start_pos, end_pos, cn_type)
                                DO UPDATE SET chrom = EXCLUDED.chrom
                                RETURNING variant_id;
                                """, 
                                [chrom, start_pos, end_pos, cn_type],True)
                        variant_cache[key]=variant_id
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
                self.batch_upsert(
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

                self.batch_upsert(
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

    def count_exons(self, region, total_exons=None):
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
            end_exon = total_exons
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
                

    def safe_json_value(self, val):
        """
        處理存進josn中的NaN 
        2026/3/20 update: deal with bool
        """
        if isinstance(val, dict):
            return {k: self.safe_json_value(v) for k, v in val.items()}
        if isinstance(val, list):
            return [self.safe_json_value(v) for v in val]
        if pd.isna(val):
            return None
        if isinstance(val, bool):
            return bool(val)
        if isinstance(val, (np.floating, float)):
            return float(val)
        if isinstance(val, (np.integer, int)):
            return int(val)
        return val
                
    def write_done_file(self, output_path):
        databasefile = os.path.join(output_path, "database.txt")
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with open(databasefile, "a") as f:
            f.write(f"Import to database completed on {now}\n")

    def variant_acmg_check(self):
        df_variant = self.sqlquery("""
                SELECT v.variant_id, v.chrom, v.pos, v.ref_base, v.alt_base, v.rsid
                FROM variant v
                LEFT JOIN variant_acmg a
                    ON v.variant_id = a.variant_id
                WHERE a.variant_id IS NULL;
        """)
        df_variant = self.add_intervar_pubmedid(df_variant)
        self.import_variant_acmg(df_variant)

    def add_intervar_pubmedid(self, df):
        exclude_keys = {
            "Build", "Chromosome", "Position",
            "Ref_allele", "Alt_allele", "Gene", "Intervar"
        }
        intervar_list = []
        acmg_rule_list = []
        pubmed_results = []
        for i in range(len(df)):
            chr_ = str(df.iloc[i]["chrom"])
            pos  = int(df.iloc[i]["pos"])
            ref  = df.iloc[i]["ref_base"]
            alt  = df.iloc[i]["alt_base"]
            rsid = df.iloc[i]["rsid"]
            rsid = df.iloc[i]["rsid"]
            if pd.isna(rsid):
                pmids = []
            else:
                pmids = self.query_rsid(str(rsid))
            try:
                res = self.ACMG_predict_wintervar(
                    chr_, pos, ref, alt, "hg38"
                )
                intervar_list.append(res.get("Intervar"))
                acmg_rules = {
                    k: v for k, v in res.items()
                    if k not in exclude_keys
                }
                acmg_str = ";".join(
                    f"{k}:{v}" for k, v in acmg_rules.items()
                )
                acmg_rule_list.append(acmg_str)
                if pmids:
                    pubmed_results.append(",".join(map(str, pmids)))
                else:
                    pubmed_results.append(None)
                time.sleep(0.2)
            except Exception:
                intervar_list.append(None)
                acmg_rule_list.append(None)
                pubmed_results.append(None)

        df["Intervar"] = intervar_list
        df["Intervar_ACMG"] = acmg_rule_list
        df["LitVar_PubMed_IDs"] = pubmed_results
        return df

    def ACMG_predict_wintervar(self, chr, pos, ref, alt, hg):
        chr = str(chr).replace("chr", "").replace("CHR", "").replace("Chr", "")
        url3 = (
            "http://wintervar.wglab.org/api_new.php"
            "?queryType=position"
            f"&chr={chr}"
            f"&pos={pos}"
            f"&ref={ref}"
            f"&alt={alt}"
            f"&build={hg}"
        )
        with req.urlopen(url3, timeout=10) as res3:
            data3 = json.load(res3)
            return data3

    def query_rsid(self, rsid):
        try:
            r = requests.get(f'https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/get/litvar@{rsid}%23%23/publications')
            if r.status_code == 200:
                return r.json().get('pmids',[])
            else:
                return []
        except requests.exceptions.RequestException as e:
            print(f"[LitVar skipped] rsid={rsid} ({e})")
            return []

def main():
    output_path = sys.argv[1]
    db_info = {
        'host': "172.17.0.1",
        'port': 5434,
        'database': "dbuser",
        'user': "dbuser",
        'password': "12345678"
    }
    importer = ImportToPostgres(db_info)
    try:
        importer.main_flow(output_path)
    finally:
        importer.close()

if __name__ == "__main__":
    main()
