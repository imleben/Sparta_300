# alpha_meta_combined.tsv 파일에서 subject별로 3개의 행(즉, 3반복)을 갖는 subject만 남기기
# Friedman 검정을 위한 전처리 작업

import pandas as pd

# 1. 파일 경로 설정(macOS)
input_file = "/Users/inseonghwang/onedrive/Sparta_300/05_diversity_20FNS/alpha_meta_combined/alpha_meta_combined.tsv"
output_file = "/Users/inseonghwang/onedrive/Sparta_300/05_diversity_20FNS/alpha_meta_combined/alpha_meta_combined_friedman.tsv"

# PC에서 실행할 때는 아래 경로로 변경
# input_file = "/mnt/d/onedrive/Sparta_300/05_diversity_20FNS/alpha_meta_combined/alpha_meta_combined.tsv"
# output_file = "/mnt/d/onedrive/Sparta_300/05_diversity_20FNS/alpha_meta_combined/alpha_meta_combined_friedman.tsv"

# 2. 파일 읽기
df = pd.read_csv(input_file, sep='\t')

# 3. subject별로 그룹화하여, 3개의 행(즉, 3반복)을 갖는 subject만 남기기
filtered_df = df.groupby('subject', group_keys=False).filter(lambda group: len(group) == 3)

# 4. 결과 저장
filtered_df.to_csv(output_file, sep='\t', index=False)
print(f"✅ 필터링된 파일이 저장되었습니다: {output_file}")
