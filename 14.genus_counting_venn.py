#!/usr/bin/env python3
import csv

def dedup(seq):
    seen = set()
    deduped = []
    for item in seq:
        if item not in seen:
            seen.add(item)
            deduped.append(item)
    return deduped

# 파일 경로 설정
input_file = "/mnt/d/onedrive/sparta_300/07_collapse/genus_20FNS-flt2-RF.tsv"
output_file = "/mnt/d/onedrive/sparta_300/07_collapse/genus_20FNS-flt2-RF-vdg.tsv"

# 파일 전체 라인 읽기
with open(input_file, 'r', newline='') as infile:
    lines = infile.readlines()

# 첫 줄이 "# "로 시작하는 설명문이면(단, "#OTU ID"는 taxonomy 컬럼이므로 유지) 제거
if lines:
    first_line = lines[0]
    if first_line.startswith("# ") and not first_line.startswith("#OTU ID"):
        lines = lines[1:]

# TSV 파일을 csv.reader로 읽기
reader = csv.reader(lines, delimiter='\t')
rows = list(reader)

if not rows:
    raise ValueError("파일에 데이터가 없습니다.")

header = rows[0]

# DQ_, MW_, TC_로 시작하는 샘플 컬럼의 인덱스 찾기
dq_indices = [i for i, col in enumerate(header) if col.startswith("DQ_")]
mw_indices = [i for i, col in enumerate(header) if col.startswith("MW_")]
tc_indices = [i for i, col in enumerate(header) if col.startswith("TC_")]

# taxonomy 정보는 첫 번째 컬럼(헤더가 "#OTU ID")에 있다고 가정
tax_col_index = 0

# 결과를 담을 리스트 초기화
dq_list = []
mw_list = []
tc_list = []

# 데이터 행 처리 (헤더 이후)
for row in rows[1:]:
    if len(row) <= tax_col_index:
        continue
    taxonomy_str = row[tax_col_index]
    # 세미콜론으로 구분
    tax_items = taxonomy_str.split(';')
    # 6번째 아이템이 없으면 넘어감
    if len(tax_items) < 6:
        continue
    # 6번째 아이템 추출 및 공백 제거
    genus_field = tax_items[5].strip()
    # 6번째 아이템이 "g__"로 시작하지 않거나, "g__" 또는 "__"만 있는 경우는 제외
    if not genus_field.startswith("g__") or genus_field in ["g__", "__"]:
        continue
    # g__ 접두어 제거하여 genus 이름 추출
    genus_name = genus_field.replace("g__", "", 1).strip()
    
    # 각 샘플 컬럼에 대해 RA 수치가 0이 아닌 경우 genus 이름을 해당 리스트에 추가
    for i in dq_indices:
        try:
            if float(row[i]) != 0:
                dq_list.append(genus_name)
        except ValueError:
            pass
    for i in mw_indices:
        try:
            if float(row[i]) != 0:
                mw_list.append(genus_name)
        except ValueError:
            pass
    for i in tc_indices:
        try:
            if float(row[i]) != 0:
                tc_list.append(genus_name)
        except ValueError:
            pass

# 각 컬럼별로 중복되는 문자열 제거 (순서 유지)
dq_list = dedup(dq_list)
mw_list = dedup(mw_list)
tc_list = dedup(tc_list)

# 세 개의 리스트 길이가 다를 수 있으므로 최대 길이에 맞춰 행 생성
max_len = max(len(dq_list), len(mw_list), len(tc_list))
new_table = []
# 새 헤더 작성
new_table.append(["DQ", "MW", "TC"])

# 각 인덱스별로 값 채우기 (없으면 빈 문자열)
for i in range(max_len):
    dq_val = dq_list[i] if i < len(dq_list) else ""
    mw_val = mw_list[i] if i < len(mw_list) else ""
    tc_val = tc_list[i] if i < len(tc_list) else ""
    new_table.append([dq_val, mw_val, tc_val])

# 결과를 TSV 파일로 저장
with open(output_file, 'w', newline='') as outfile:
    writer = csv.writer(outfile, delimiter='\t')
    writer.writerows(new_table)

print("변환 완료. 결과는 다음 파일에 저장되었습니다:")
print(output_file)
