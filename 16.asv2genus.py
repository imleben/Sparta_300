#  ASV ID와 Genus 정보가 포함된 메타데이터 파일을 읽어서 FASTA 파일의 헤더를 Genus로 변경하는 스크립트
#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def read_mapping_file(mapping_file):
    """메타데이터 파일에서 ASV ID와 Genus 매핑 정보 읽기"""
    mapping_df = pd.read_csv(mapping_file, sep='\t')
    # ASV ID와 Genus 컬럼만 선택하여 딕셔너리로 변환
    id_to_genus = dict(zip(mapping_df['FeatureID'], mapping_df['Genus']))
    return id_to_genus

def rename_fasta_headers(input_fasta, output_fasta, id_to_genus):
    """FASTA 파일의 헤더를 Genus 명으로 변경하고 시퀀스를 한 줄로 출력"""
    
    with open(output_fasta, 'w') as out_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            original_id = record.id.strip()
            genus = id_to_genus.get(original_id, original_id)
            
            # 헤더와 시퀀스를 직접 작성 
            out_handle.write(f">{genus}\n")
            out_handle.write(f"{str(record.seq)}\n")

def main():
    # 파일 경로 설정
    mapping_file = "/mnt/d/OneDrive/Sparta_300/10_phylo/taxa_20FNS-flt2/selected_ASVs.txt"
    input_fasta = "/mnt/d/OneDrive/Sparta_300/10_phylo/rep-seqs_20FNS-flt2/genus-repseq_exported/dna-sequences.fasta"
    output_fasta = "/mnt/d/OneDrive/Sparta_300/10_phylo/rep-seqs_20FNS-flt2/genus-repseq_exported/dna-sequences_genus_mapped.fasta"
    # ASV ID -> Genus 매핑 읽기
    print("매핑 파일 읽는 중...")
    id_to_genus = read_mapping_file(mapping_file)
    print(f"총 {len(id_to_genus)} 개의 매핑 정보를 읽었습니다.")
    
    # FASTA 파일 처리
    print("FASTA 파일 처리 중...")
    rename_fasta_headers(input_fasta, output_fasta, id_to_genus)
    print(f"처리 완료. 결과가 {output_fasta}에 저장되었습니다.")

if __name__ == "__main__":
    main()
