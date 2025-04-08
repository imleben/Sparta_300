#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import re
import os
import sys

def check_files_exist(taxonomy_file, feature_table_file):
    tax_exists = os.path.exists(taxonomy_file)
    feat_exists = os.path.exists(feature_table_file)
    
    print(f"taxonomy 파일 존재: {tax_exists}")
    print(f"feature-table 파일 존재: {feat_exists}")
    
    return tax_exists and feat_exists

def read_taxonomy_file(taxonomy_file):
    print(f"\n{taxonomy_file} 파일 읽는 중...")
    
    # 다양한 구분자 시도
    for sep in ['\t', ',']:
        try:
            taxonomy_df = pd.read_csv(taxonomy_file, sep=sep)
            if len(taxonomy_df.columns) > 1:
                break
        except:
            continue
    
    print(f"▶ Taxonomy 파일 형태: {taxonomy_df.shape}")
    print(f"  Taxonomy 열 이름: {list(taxonomy_df.columns)}")
    print("  Taxonomy 샘플 데이터:")
    print(taxonomy_df.head(2))
    
    # Feature ID 열 이름 확인
    feature_id_col = None
    for col in taxonomy_df.columns:
        if 'id' in col.lower() or 'feature' in col.lower() or 'otu' in col.lower():
            feature_id_col = col
            break
    
    if not feature_id_col:
        feature_id_col = taxonomy_df.columns[0]  # 첫 번째 열을 ID로 가정
    
    print(f"▶ Feature ID 형식: {taxonomy_df[feature_id_col].iloc[0]}")
    
    return taxonomy_df, feature_id_col

def read_feature_table(feature_table_file):
    print(f"\n{feature_table_file} 파일 읽는 중...")
    
    # 파일 먼저 검사
    with open(feature_table_file, 'r') as file:
        first_lines = [file.readline() for _ in range(5)]
    
    print("▶ Feature Table 파일 첫 몇 줄:")
    for line in first_lines[:2]:
        print(f"  {line.strip()}")
    
    # 헤더 라인 식별
    header_line = 0
    for i, line in enumerate(first_lines):
        if '#OTU ID' in line or 'Feature ID' in line:
            header_line = i
            break
    
    print(f"▶ 헤더 라인 위치: 라인 {header_line}")
    
    # 파일 형식에 따라 파싱 방법 결정
    try:
        # 방법 1: pandas로 직접 읽기
        feature_table = pd.read_csv(feature_table_file, sep='\t', skiprows=header_line)
    except Exception as e:
        print(f"표준 파싱 실패: {e}")
        # 방법 2: 수동 파싱
        with open(feature_table_file, 'r') as file:
            lines = file.readlines()[header_line:]
        
        header = lines[0].strip().split('\t')
        data = []
        for line in lines[1:]:
            values = line.strip().split('\t')
            if len(values) > 1:
                row = {header[i]: values[i] for i in range(min(len(header), len(values)))}
                data.append(row)
        
        feature_table = pd.DataFrame(data)
    
    print(f"▶ Feature Table 형태: {feature_table.shape}")
    
    # 첫 번째 열 이름 식별
    feature_id_col = feature_table.columns[0]
    print(f"  Feature Table 첫 열 이름: {feature_id_col}")
    print("  Feature Table 샘플 데이터:")
    print(feature_table.head(1))
    
    return feature_table, feature_id_col

def process_feature_table(feature_table, feature_id_col):
    # 빈도 데이터프레임 생성
    # 첫 번째 열을 Feature ID로, 나머지 열의 합계를 Total Frequency로
    try:
        # 숫자 열만 선택하여 합계 계산
        numeric_cols = feature_table.columns[1:]
        feature_table[numeric_cols] = feature_table[numeric_cols].apply(pd.to_numeric, errors='coerce')
        
        frequency_df = pd.DataFrame({
            'FeatureID': feature_table[feature_id_col],
            'TotalFrequency': feature_table[numeric_cols].sum(axis=1)
        })
    except Exception as e:
        print(f"빈도 계산 오류: {e}")
        # 대안적 방법
        print("대안적 빈도 계산 시도 중...")
        
        # 첫 번째 열을 제외한 모든 열을 문자열 -> 숫자로 변환
        sample_cols = feature_table.columns[1:]
        frequency_data = []
        
        for _, row in feature_table.iterrows():
            feature_id = row[feature_id_col]
            frequencies = []
            
            for col in sample_cols:
                try:
                    val = float(row[col])
                    frequencies.append(val)
                except:
                    frequencies.append(0)
            
            total_freq = sum(frequencies)
            frequency_data.append({'FeatureID': feature_id, 'TotalFrequency': total_freq})
        
        frequency_df = pd.DataFrame(frequency_data)
    
    print(f"▶ 빈도 데이터프레임 형태: {frequency_df.shape}")
    print("  빈도 데이터프레임 샘플:")
    print(frequency_df.head(2))
    
    # Feature ID 샘플 확인
    feature_ids = frequency_df['FeatureID'].head(5).tolist()
    print(f"▶ Feature IDs 샘플: {feature_ids}")
    
    return frequency_df

def merge_dataframes(taxonomy_df, frequency_df, tax_id_col, freq_id_col):
    print("\n데이터프레임 병합 중...")
    
    print(f"▶ Taxonomy ID 형식: {taxonomy_df[tax_id_col].iloc[0]}")
    print(f"  Frequency ID 형식: {frequency_df['FeatureID'].iloc[0]}")
    
    # ID 열 이름 통일
    taxonomy_df = taxonomy_df.rename(columns={tax_id_col: 'FeatureID'})
    print("▶ ID 열 이름 통일 후:")
    print(f"  Taxonomy FeatureID 열: FeatureID")
    print(f"  Frequency FeatureID 열: FeatureID")
    
    # ID 형식이 다를 경우 변환 시도
    taxonomy_ids = set(taxonomy_df['FeatureID'])
    frequency_ids = set(frequency_df['FeatureID'])
    common_ids = taxonomy_ids.intersection(frequency_ids)
    
    print(f"  Taxonomy ID 샘플: {list(taxonomy_df['FeatureID'].head(5))}")
    print(f"  Frequency ID 샘플: {list(frequency_df['FeatureID'].head(5))}")
    print(f"  공통 ID 개수: {len(common_ids)}")
    
    if len(common_ids) == 0:
        print("▶ ID 형식이 다릅니다. 변환 시도 중...")
        
        # 가능한 변환 시도
        # 예: DNA 서열에서 해시로, 또는 다른 일반적인 변환 패턴 적용
        
        # 1. 해시 ID -> 서열 번호 변환 시도
        try:
            # taxonomy는 해시 ID, frequency는 번호 ID일 경우
            map_dict = {f"ASV{i+1}": id for i, id in enumerate(taxonomy_df['FeatureID'])}
            frequency_df['HashID'] = frequency_df['FeatureID'].map(map_dict)
            merged_df = pd.merge(taxonomy_df, frequency_df, left_on='FeatureID', right_on='HashID')
            if len(merged_df) > 0:
                print("▶ 해시 ID -> 서열 번호 변환 성공")
                return merged_df
        except:
            pass
        
        # 2. 번호 ID -> 해시 ID 변환 시도
        try:
            # taxonomy는 번호 ID, frequency는 해시 ID일 경우
            map_dict = {id: f"ASV{i+1}" for i, id in enumerate(frequency_df['FeatureID'])}
            taxonomy_df['NumID'] = taxonomy_df['FeatureID'].map(map_dict)
            merged_df = pd.merge(taxonomy_df, frequency_df, left_on='NumID', right_on='FeatureID')
            if len(merged_df) > 0:
                print("▶ 번호 ID -> 해시 ID 변환 성공")
                return merged_df
        except:
            pass
            
        print("▶ 자동 변환 실패. 수동으로 ID 매핑 필요")
        return pd.DataFrame()
    else:
        # 직접 병합
        merged_df = pd.merge(taxonomy_df, frequency_df, on='FeatureID')
        
    print(f"▶ 병합된 데이터프레임 형태: {merged_df.shape}")
    
    if len(merged_df) > 0:
        print("  병합된 데이터프레임 샘플:")
        print(merged_df.head(2))
    else:
        print("병합 결과 없음")
        print("▶ 병합 결과가 없습니다. Feature ID 형식 불일치 문제를 확인하세요.")
    
    return merged_df

def extract_genus(merged_df, taxon_col='Taxon'):
    print("\nGenus 추출 중...")
    
    if len(merged_df) == 0:
        return pd.DataFrame()
    
    # Taxon 열에서 Genus 추출
    def extract_genus_from_taxon(taxon):
        try:
            # QIIME2 형식 (d__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus)
            if ';' in taxon:
                genus_match = re.search(r'g__([^;]+)', taxon)
                if genus_match:
                    return genus_match.group(1).strip()
            
            # 또는 다른 형식 (Bacteria;Firmicutes;Bacilli;Lactobacillales;Streptococcaceae;Streptococcus)
            parts = taxon.split(';')
            if len(parts) >= 6:
                return parts[5].strip()
            
            return "Unknown"
        except:
            return "Unknown"
    
    genus_df = pd.DataFrame({
        'Taxon': merged_df[taxon_col],
        'Genus': merged_df[taxon_col].apply(extract_genus_from_taxon)
    })
    
    print("추출된 Genus 샘플:")
    print(genus_df.head(5))
    
    return genus_df

def select_representative_asvs(merged_df, genus_df):
    print("\n각 Genus별 최대 빈도 Feature ID 선택 중...")
    
    if len(merged_df) == 0 or len(genus_df) == 0:
        return pd.DataFrame()
    
    # Genus 열을 merged_df에 추가
    merged_df['Genus'] = genus_df['Genus']
    
    # 각 Genus별로 가장 빈도가 높은 ASV 선택
    selected_asvs = merged_df.loc[merged_df.groupby('Genus')['TotalFrequency'].idxmax()]
    
    print(f"선택된 데이터프레임 형태: {selected_asvs.shape}")
    
    if len(selected_asvs) > 0:
        print("선택된 데이터프레임 샘플:")
        print(selected_asvs.head(2))
    
    return selected_asvs

def save_results(selected_asvs, output_file):
    if len(selected_asvs) == 0:
        print("\n선택된 Genus가 없습니다. 데이터를 확인해 주세요.")
        return
    
    # 출력 디렉토리 확인
    output_dir = os.path.dirname(output_file)
    if not os.path.exists(output_dir):
        response = input(f"\n저장 경로 {output_dir}가 존재하지 않습니다. 생성하시겠습니까? (y/n): ")
        if response.lower() == 'y':
            try:
                os.makedirs(output_dir)
                print(f"디렉토리 생성됨: {output_dir}")
            except Exception as e:
                print(f"디렉토리 생성 실패: {e}")
                return
        else:
            print("작업이 취소되었습니다.")
            return
    
    # 결과 저장
    selected_asvs.to_csv(output_file, sep='\t', index=False)
    print(f"\n{output_file} 파일에 결과가 저장되었습니다.")
    print(f"총 {len(selected_asvs)} 개의 Genus에서 대표 ASV를 선택했습니다.")


def main():
    # 절대 경로 사용
    taxonomy_file = "/mnt/d/OneDrive/Sparta_300/10_phylo/taxa_20FNS/taxonomy.tsv"
    feature_table_file = "/mnt/d/OneDrive/Sparta_300/10_phylo/table_20FNS-flt3/genus_20FNS_feature_table.tsv"
    output_file = "/mnt/d/OneDrive/Sparta_300/10_phylo/taxa_20FNS-flt3/selected_ASVs.txt"

    
    # 파일 존재 확인
    if not check_files_exist(taxonomy_file, feature_table_file):
        print("파일이 존재하지 않습니다. 경로를 확인해주세요.")
        return
    
    # Taxonomy 파일 읽기
    taxonomy_df, tax_id_col = read_taxonomy_file(taxonomy_file)
    
    # Feature-table 파일 읽기
    feature_table, feature_id_col = read_feature_table(feature_table_file)
    
    # Feature 빈도 처리
    frequency_df = process_feature_table(feature_table, feature_id_col)
    
    # 데이터프레임 병합
    merged_df = merge_dataframes(taxonomy_df, frequency_df, tax_id_col, 'FeatureID')
    
    # Genus 추출
    genus_df = extract_genus(merged_df)
    
    # 대표 ASV 선택
    selected_asvs = select_representative_asvs(merged_df, genus_df)
    
    # 결과 저장
    save_results(selected_asvs, output_file)

if __name__ == "__main__":
    main()
