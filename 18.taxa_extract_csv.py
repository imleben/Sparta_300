# Description: OTU 테이블에서 분류학적 정보를 추출하여 CSV 파일로 저장하고, 원본 파일을 CSV로 변환하며 첫 열을 추출된 OTU로 대체한 후 전치된 파일을 생성하는 스크립트. 전치는 필요할 경우에만 사용.

import pandas as pd
import re
import os

# 기본 디렉토리 경로 한 번만 설정
base_dir = '/Users/inseonghwang/OneDrive/Sparta_300/07_collapse/biom_genus_flt1'
base_filename = 'genus_20FNS-flt1'

# 파일 경로 자동 생성
input_file = os.path.join(base_dir, f"{base_filename}.tsv")
output_taxonomy_file = os.path.join(base_dir, f"{base_filename}-taxa-only.csv")
output_converted_file = os.path.join(base_dir, f"{base_filename}.csv")
output_transposed_file = os.path.join(base_dir, f"{base_filename}_transposed.csv")

# 1. 먼저 taxonomy 정보 추출하여 저장
def extract_taxonomy():
    print(f"Step 1: Extracting taxonomy information from {input_file}")
    # 첫 번째 컬럼만 추출하여 임시 저장
    with open(input_file, 'r') as file:
        lines = file.readlines()
        
        # 첫 번째 컬럼만 추출
        otu_ids = []
        for line in lines:
            if line.startswith('#OTU ID') or line.startswith('d__'):
                parts = line.strip().split('\t')
                otu_ids.append(parts[0])

    # 분류학적 정보 추출 함수
    def extract_taxonomy_levels(otu_id):
        # 정규식으로 각 분류 수준 추출
        domain_match = re.search(r'd__([^;]+)', otu_id)
        phylum_match = re.search(r'p__([^;]+)', otu_id)
        class_match = re.search(r'c__([^;]+)', otu_id)
        order_match = re.search(r'o__([^;]+)', otu_id)
        family_match = re.search(r'f__([^;]+)', otu_id)
        genus_match = re.search(r'g__([^;]+)', otu_id)
        
        # 매치된 결과에서 그룹 추출 또는 빈 문자열 처리
        domain = domain_match.group(1) if domain_match else ""
        phylum = phylum_match.group(1) if phylum_match else ""
        class_ = class_match.group(1) if class_match else ""
        order = order_match.group(1) if order_match else ""
        family = family_match.group(1) if family_match else ""
        
        # Genus 정보 추출 (숫자나 공백 이전까지의 문자열만)
        genus = ""
        if genus_match:
            genus_full = genus_match.group(1)
            genus = genus_full.split()[0] if ' ' in genus_full else genus_full
        
        return {
            "Domain": domain,
            "Phylum": phylum,
            "Class": class_,
            "Order": order,
            "Family": family,
            "Genus": genus
        }

    # 데이터 추출
    data = []
    for otu_id in otu_ids[1:]:  # 첫 번째 줄(헤더)은 건너뜀
        taxonomy_info = extract_taxonomy_levels(otu_id)
        
        # g__가 없는 경우에도 행을 유지하며, 적절한 대체값 사용
        if taxonomy_info['Genus'] == "":
            taxonomy_level = "Unknown"
            # 가장 낮은 수준의 분류학적 정보를 찾아 대체값 결정
            if taxonomy_info['Family']:
                taxonomy_level = f"f__{taxonomy_info['Family']}"
            elif taxonomy_info['Order']:
                taxonomy_level = f"o__{taxonomy_info['Order']}"
            elif taxonomy_info['Class']:
                taxonomy_level = f"c__{taxonomy_info['Class']}"
            elif taxonomy_info['Phylum']:
                taxonomy_level = f"p__{taxonomy_info['Phylum']}"
            elif taxonomy_info['Domain']:
                taxonomy_level = f"d__{taxonomy_info['Domain']}"
            
            # #TAXONOMY 열에는 "Unclassified_[가장 낮은 분류학적 수준]"으로 표시
            data.append({
                "#TAXONOMY": f"Unclassified_{taxonomy_level}",
                **taxonomy_info
            })
        else:
            data.append({
                "#TAXONOMY": f"g__{taxonomy_info['Genus']}",
                **taxonomy_info
            })

    # DataFrame 생성 및 저장
    df = pd.DataFrame(data)
    df.to_csv(output_taxonomy_file, index=False)

    print(f"Taxonomy data extracted and saved to {output_taxonomy_file}")
    print(f"Total records: {len(data)}")
    
    # 각 분류 수준별 레코드 수 확인
    print("\nTaxonomy level counts:")
    for level in ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus']:
        valid_count = df[df[level] != ''].shape[0]
        print(f"{level}: {valid_count} valid entries ({valid_count/len(data)*100:.1f}%)")
    
    return df

# 2. 원본 TSV 파일을 CSV로 변환하고, 첫 컬럼을 추출된 taxonomy의 첫 컬럼으로 대체
def convert_to_csv(taxonomy_df):
    print(f"\nStep 2: Converting original TSV file to CSV and replacing first column")
    
    # 원본 파일을 TSV로 읽음
    with open(input_file, 'r') as file:
        lines = file.readlines()
    
    # 헤더와 데이터 행 분리
    header = lines[0].strip().split('\t')
    data_rows = [line.strip().split('\t') for line in lines[1:]]
    
    # 추출된 taxonomy의 첫 컬럼 가져오기 ("#TAXONOMY" 열)
    taxonomy_column = taxonomy_df['#TAXONOMY'].tolist()
    
    # 새 CSV 파일 생성
    with open(output_converted_file, 'w') as out_file:
        # 헤더의 첫 번째 열을 "#NAME"으로 변경
        header[0] = "#NAME"
        out_file.write(','.join(header) + '\n')
        
        # 각 행에 대해 첫 번째 열을 taxonomy 정보로 대체
        for i, row in enumerate(data_rows):
            if i < len(taxonomy_column):  # 인덱스 범위 체크
                row[0] = taxonomy_column[i]
            
            # CSV 형식으로 쓰기 (따옴표 처리 포함)
            csv_row = []
            for item in row:
                # 쉼표가 포함된 항목은 따옴표로 감싸기
                if ',' in item:
                    csv_row.append(f'"{item}"')
                else:
                    csv_row.append(item)
            
            out_file.write(','.join(csv_row) + '\n')
    
    print(f"Original file converted to CSV with modified first column: {output_converted_file}")
    
    # 변환된 파일 미리보기
    preview_df = pd.read_csv(output_converted_file, nrows=5)
    print("\nPreview of the converted CSV file:")
    print(preview_df.head())
    
    # 3. 전치된 파일 생성
    print("\nStep 3: Creating transposed version of the CSV file")
    transposed_df = pd.read_csv(output_converted_file, index_col=0).T
    
    # 전치된 DataFrame에서 인덱스 이름 설정 (이전 컬럼 이름이 인덱스가 됨)
    transposed_df.index.name = "#NAME"
    
    # 전치된 파일 저장
    transposed_df.to_csv(output_transposed_file)
    
    print(f"Transposed file created and saved to {output_transposed_file}")
    print("\nPreview of the transposed CSV file:")
    print(transposed_df.head())
    
    return transposed_df

# 메인 실행 부분
def main():
    # 파일 존재 확인
    if not os.path.exists(input_file):
        print(f"Error: Input file not found at {input_file}")
        return
    
    # 1. Taxonomy 추출 및 저장
    taxonomy_df = extract_taxonomy()
    
    # 2. TSV를 CSV로 변환하면서 첫 열 대체 및 전치 파일 생성
    transposed_df = convert_to_csv(taxonomy_df)
    
    print("\nProcess completed successfully. Three files created:")
    print(f"1. Taxonomy file: {output_taxonomy_file}")
    print(f"2. Converted CSV: {output_converted_file}")
    print(f"3. Transposed CSV: {output_transposed_file}")

if __name__ == "__main__":
    main()