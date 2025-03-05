import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# ============================================
# 1. 경로 설정 
# ============================================
input_path = "/Users/inseonghwang/OneDrive/Sparta_300/08_heatmap_20FNS/genus_20FNS-flt2.tsv"
output_folder = "/Users/inseonghwang/OneDrive/Sparta_300/08_heatmap_20FNS"
output_filename = "genus_heatmap_20FNS-flt2_top30.pdf"
output_path = f"{output_folder}/{output_filename}"

# ============================================
# 2. 데이터 읽기 (relative abundance 변환 전 파일 이용)
# ============================================
abundances = pd.read_table(input_path, skiprows=1, index_col=0)

# 데이터 복사
abund_to_plot = abundances.copy()

# Genus 이름만 사용하기 (세미콜론으로 구분된 문자열의 6번째 항목 사용)
abund_to_plot.index = abund_to_plot.index.str.split(";").str[5]

# "unclassified" 또는 잘못된 index 삭제
abund_to_plot = abund_to_plot[~abund_to_plot.index.isin(["g__", "__"])]

# 각 genus의 전체 abundance 합을 계산하고, 많은 순으로 정렬하여 상위 30개 선택
abund_to_plot['total_abundance'] = abund_to_plot.sum(axis=1)
abund_to_plot = abund_to_plot.sort_values(by='total_abundance', ascending=False)
abund_to_plot = abund_to_plot.head(30)
abund_to_plot = abund_to_plot.drop(columns=['total_abundance'])

# ============================================
# 3. Centered Log-Ratio (CLR) 변환 함수 및 적용
# ============================================
def clr_transform(xs):
    logged = np.log(xs + 0.5)
    return logged - logged.mean()

transformed = abund_to_plot.apply(clr_transform, axis=0)

# ============================================
# 4. 클러스터 맵 그리기
# ============================================
cluster_grid = sns.clustermap(transformed.T, cmap="magma", xticklabels=True, figsize=(8, 15))

# ============================================
# 5. 그래프를 파일로 저장
# ============================================
cluster_grid.savefig(output_path)

print(f"Heatmap saved as '{output_filename}' in {output_folder}")
