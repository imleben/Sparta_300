############################################################
# HOMD 분류체계 파일을 Silva 형식으로 변환하고 공통 식별자 생성
############################################################

import pandas as pd
import re

def convert_homd_to_silva(input_file, output_file):
    """
    HOMD 분류체계 파일을 Silva 형식으로 변환하고 공통 식별자 생성
    
    Args:
        input_file (str): HOMD taxonomy 파일 경로
        output_file (str): 출력할 Silva 형식 파일 경로
    """
    # 파일 읽기 (첫 두 줄은 건너뛰기)
    with open(input_file, 'r', encoding='utf-8') as f:
        # 첫 두 줄 건너뛰기
        next(f)  # 파일명 줄
        header = next(f).strip().split('\t')  # 헤더 줄
        
        # 데이터 읽기
        data = []
        for line in f:
            values = line.strip().split('\t')
            if len(values) >= 9:  # 최소한 Species 열까지 존재하는지 확인
                row_data = dict(zip(header, values))
                data.append(row_data)
    
    # 데이터프레임 생성
    df = pd.DataFrame(data)
    
    # Silva 형식으로 변환
    silva_format = []
    
    for _, row in df.iterrows():
        try:
            # HMT_ID에 'HMT-' 접두사 추가하여 공통 식별자 생성
            hmt_id = row['HMT_ID']
            common_id = f"HMT-{hmt_id}"
            
            # Species 이름 처리
            species = row['Species']
            genus = row['Genus']
            
            # silva에서는 보통 genus_species 형식 사용
            if species.startswith('sp._'):  # sp._HMT_XXX 형식인 경우
                species_format = species.replace('.', '')  # 점 제거
            else:
                species_format = f"{genus}_{species}"
            
            # 분류체계 문자열 생성
            taxonomy = [
                f"d__{row['Domain']}",
                f"p__{row['Phylum']}",
                f"c__{row['Class']}",
                f"o__{row['Order']}",
                f"f__{row['Family']}",
                f"g__{genus}",
                f"s__{species_format}"
            ]
            
            silva_format.append({
                'Feature ID': common_id,  # 'HMT-XXX' 형식의 공통 식별자 사용
                'Taxon': '; '.join(taxonomy)
            })
        except KeyError as e:
            print(f"키 오류 발생: {e}. 행 건너뛰기: {row}")
            continue
    
    # 결과 데이터프레임 생성 및 저장
    result_df = pd.DataFrame(silva_format)
    result_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"변환 완료: {len(silva_format)}개 항목이 {output_file}에 저장되었습니다.")
    print(f"HMT-XXX 형식의 공통 식별자가 Feature ID 열에 추가되었습니다.")

# 스크립트 실행 예시
if __name__ == "__main__":
    input_file = "HOMD_taxon_table2025-03-13_1741906836.txt"
    output_file = "homd2silva_taxon.tsv"
    
    convert_homd_to_silva(input_file, output_file)
