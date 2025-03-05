#!/usr/bin/env python3
"""
기능 요약:
  1. 입력 파일에서 주석(#) 행 제거
     - "#OTU ID"로 시작하지 않는 모든 # 주석 행 제거
  2. 헤더(`#OTU ID`)를 그대로 둔 상태로, "class" 행 생성 후 맨 위(0번 인덱스)에 삽입
     - "DQ", "MW", "TC"로 시작하는 열만 접두어(prefix) 추출
     - 그렇지 않은 열은 공백("") 추가
  3. 원래 헤더 행(#OTU ID 행)을 복사하여 두 행으로 만듦
     - 첫 번째 헤더 행의 첫 셀은 "subclass"로 변경
     - 복사된(두 번째) 헤더 행의 첫 셀은 "subject"로 변경
  4. 나머지 데이터는 그대로 유지
  5. 전치 기능(기존 코드)은 **주석 처리**해 필요 시 주석 해제 후 사용
"""

import csv
import sys
import os
import pandas as pd

# --- (1) 수정하기 쉽게 코드 상단에 경로를 배치 ---
# INPUT_FILE = "/mnt/d/onedrive/sparta_300/taxa_20FNS/exported_species_RF/species_20FNS_fltrd_RF.tsv"
# OUTPUT_FILE = "/mnt/d/onedrive/sparta_300/taxa_20FNS/LEfSe/input_table.tsv"

INPUT_FILE = "/Users/inseonghwang/onedrive/sparta_300/taxa_20FNS/exported_species_RF/species_20FNS_fltrd_RF.tsv"
OUTPUT_FILE = "/Users/inseonghwang/onedrive/sparta_300/taxa_20FNS/LEfSe/input_table.tsv"

def insert_class_subclass_sampleid_rows(input_path, output_path):
    # 1) 파일 전체 읽기
    with open(input_path, "r", newline="") as infile:
        reader = csv.reader(infile, delimiter="\t")
        rows = list(reader)

    # --- (1-A) 주석(#) 행 제거 ---
    # '#OTU ID'가 들어있는 행은 헤더로 남기고,
    # 그 외 '#'로 시작하는 행은 모두 제거
    while rows and rows[0] and rows[0][0].startswith('#') and not rows[0][0].lower().startswith('#otu id'):
        rows.pop(0)

    # 헤더 유효성 체크
    if not rows:
        sys.exit("유효한 헤더를 찾지 못했습니다. (#OTU ID로 시작하는 행이 필요합니다.)")

    # 현재 rows[0]가 원본 헤더 행이어야 함 (예: ["#OTU ID", "DQ_001", "MW_003", ...])
    original_header = rows[0]

    # --- (2) "class" 행 생성 ---
    #     첫 열은 "class", 헤더의 두 번째 열부터 "DQ"/"MW"/"TC" 접두어 추출
    class_row = ["class"]
    for cell in original_header[1:]:
        cell_str = cell.strip()
        if cell_str.startswith("DQ") or cell_str.startswith("MW") or cell_str.startswith("TC"):
            prefix = cell_str.split("_")[0]
            class_row.append(prefix)
        else:
            class_row.append("")

    # class_row를 맨 위(인덱스 0)에 삽입
    rows.insert(0, class_row)
    # 이제 rows[1]이 원본 헤더 행이 됨

    # --- (3) 원본 헤더 행을 2개로 만들어
    #         첫 번째 헤더 행 → 첫 셀 "subclass"
    #         두 번째 헤더 행(복사본) → 첫 셀 "subject"
    # 현재 rows[1]이 "#OTU ID ...",
    # 이를 수정하여 "subclass ...",
    # 그리고 복사본을 만들어 "subject ..."로 함

    # (3-A) 첫 번째 헤더 행의 첫 셀을 "subclass"로 변경
    rows[1][0] = "subclass"

    # (3-B) 두 번째 헤더 복제하여 첫 셀만 "subject"로 변경
    #       list(rows[1]) 로 하면 얕은 복사이지만 여기서는 충분
    header_copy = list(rows[1])
    header_copy[0] = "subject"

    # (3-C) 복사본을 rows[2] 위치에 삽입
    rows.insert(2, header_copy)

    # --- (4) 전치 기능은 주석 처리 ---
    df = pd.DataFrame(rows)
    """
    # 아래 주석 해제 시 전치 수행 + 필요 시 #OTU ID → subject 등 추가 로직
    df_transposed = df.T

    # 전치 후 첫 행/열 수정이 필요하면 여기에 작성

    df_transposed.to_csv(output_path, sep="\t", index=False, header=False)
    return
    """

    # (전치를 하지 않고 그대로 TSV로 저장)
    df.to_csv(output_path, sep="\t", index=False, header=False)


def main():
    try:
        insert_class_subclass_sampleid_rows(INPUT_FILE, OUTPUT_FILE)
        print(f"[완료] class/subclass/subject 행 처리 완료. (전치 기능 주석 상태)\n출력: {OUTPUT_FILE}")
    except Exception as e:
        print(f"[오류 발생] {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
