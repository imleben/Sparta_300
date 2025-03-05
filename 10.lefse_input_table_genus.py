# -*- coding: utf-8 -*-
# LEfSe 분석용 형식 변환 스크립트

input_path = "/mnt/d/onedrive/sparta_300/08_heatmap_20FNS/genus_20FNS-flt2.tsv"
output_path = "/mnt/d/onedrive/sparta_300/09_LEfSe/input_table_genus.tsv"

# 1. 입력 파일을 읽고 주석 줄('#'로 시작하는 줄)을 필터링합니다.
header_line = None
data_lines = []
with open(input_path, 'r') as infile:
    for line in infile:
        if line.startswith('#'):
            # '#OTU ID'를 포함한 헤더 라인은 유지하고 다른 주석은 무시
            if line.strip().startswith("#OTU ID"):
                header_line = line.strip()  # 헤더 행 저장 (양끝 공백 제거)
            # 그 외 '#'로 시작하는 줄은 건너뜀
            else:
                continue
        else:
            # 데이터 라인은 리스트에 추가 (개행문자 제거)
            data_lines.append(line.strip())

# 헤더가 존재하는지 확인
if header_line is None:
    raise RuntimeError("입력 파일에 '#OTU ID' 헤더 행을 찾을 수 없습니다.")

# 2. 헤더 행을 탭으로 분리하고, 각 샘플 열 이름에서 접두사 추출하여 group 행 생성
headers = header_line.split('\t')
group_row = ["group"]  # 첫 번째 열의 group 행 이름
for col_name in headers[1:]:  # 첫 열('#OTU ID') 제외한 샘플 열들의 접두사 추출
    if col_name.startswith("DQ"):
        group_row.append("DQ")
    elif col_name.startswith("MW"):
        group_row.append("MW")
    elif col_name.startswith("TC"):
        group_row.append("TC")
    else:
        # 지정된 접두사(DQ, MW, TC) 외의 열은 그대로 사용하거나 기타 처리
        group_row.append(col_name)

# 3. subject 행 생성: 원본 헤더를 사용하되 첫 번째 항목을 'subject'로 변경
subject_row = ["subject"] + headers[1:]

# 4. subgroup 행 생성 로직은 필요 시 추가

# 5. 출력 파일에 group 행, subject 행, 그리고 원본 데이터 행들을 순서대로 기록
with open(output_path, 'w', encoding='utf-8') as outfile:
    # group 행 작성
    outfile.write('\t'.join(group_row) + '\n')
    # subject 행 작성
    outfile.write('\t'.join(subject_row) + '\n')
    # 원본 데이터 행들 작성 (첫 컬럼 처리 추가)
    for line in data_lines:
        columns = line.split('\t')
        # 첫 번째 컬럼에 세미콜론(';')이 포함된 경우 처리
        if ';' in columns[0]:
            parts = columns[0].split(';')
            # 6번째 아이템(인덱스 5)이 존재하고 g__로 시작하면 그 값만 사용
            if len(parts) >= 6:
                genus = parts[5].strip()
                if genus.startswith("g__"):
                    columns[0] = genus
        outfile.write('\t'.join(columns) + '\n')

print(f"변환된 테이블이 '{output_path}' 경로에 저장되었습니다.")
