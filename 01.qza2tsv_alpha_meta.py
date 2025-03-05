import os
import subprocess
import pandas as pd
from functools import reduce

# ============================================
# 1. 경로 설정
# ============================================
input_dir = "/mnt/d/Onedrive/Sparta_300/diversity_20FNS"
metadata_file = "/mnt/d/Onedrive/Sparta_300/metadata_20FNS.tsv"
output_dir = "/mnt/d/Onedrive/Sparta_300/diversity_20FNS/alpha_meta_combined"
os.makedirs(output_dir, exist_ok=True)

# ============================================
# 2. QZA → TSV 변환 (파일명 alpha-diversity.tsv로 수정)
# ============================================
def process_qza(qza_path):
    """QZA 파일을 TSV로 변환하고 컬럼 이름 보정"""
    try:
        # 메트릭 이름 추출 (예: chao1_vector.qza → chao1)
        metric_name = os.path.basename(qza_path).replace("_vector.qza", "")
        
        # 임시 폴더에 QIIME2 export
        temp_dir = os.path.join(output_dir, f"TEMP_{metric_name}")
        subprocess.run(
            ["qiime", "tools", "export", "--input-path", qza_path, "--output-path", temp_dir],
            check=True
        )
        
        # export된 파일은 alpha-diversity.tsv 라고 가정
        temp_tsv = os.path.join(temp_dir, "alpha-diversity.tsv")
        df = pd.read_csv(temp_tsv, sep='\t')
        
        # 컬럼 이름 강제 지정: 첫 컬럼 → sample.id, 두 번째 컬럼 → metric_name
        df.columns = ['sample.id', metric_name]
        
        # 임시 폴더 삭제
        subprocess.run(["rm", "-rf", temp_dir])
        return df
    
    except Exception as e:
        print(f"[ERROR] {qza_path} 처리 실패: {str(e)}")
        return None

# ============================================
# 3. 모든 QZA 파일 처리
# ============================================
all_dfs = []
for file in os.listdir(input_dir):
    if file.endswith("_vector.qza"):
        qza_path = os.path.join(input_dir, file)
        df = process_qza(qza_path)
        if df is not None:
            all_dfs.append(df)

# ============================================
# 4. 데이터 병합 및 메타데이터 통합
# ============================================
if all_dfs:
    merged_diversity = reduce(
        lambda left, right: pd.merge(left, right, on='sample.id', how='outer'), 
        all_dfs
    )
    
    # 메타데이터 불러오기
    metadata = pd.read_csv(metadata_file, sep='\t')
    
    # 컬럼 이름에 포함된 공백 제거
    metadata.columns = metadata.columns.str.strip()
    
    # metadata의 첫 번째 컬럼 이름(어떤 이름이든 상관없음)을 'sample.id'로 강제 변경
    first_col = metadata.columns[0]
    metadata.rename(columns={first_col: 'sample.id'}, inplace=True)
    
    # merged_diversity와 metadata를 'sample.id' 컬럼을 기준으로 병합
    final_df = pd.merge(merged_diversity, metadata, on='sample.id', how='left')
    
    output_path = os.path.join(output_dir, "alpha_meta_combined.tsv")
    final_df.to_csv(output_path, sep='\t', index=False)
    print(f"\n✅ 최종 파일 저장됨: {output_path}")
else:
    print("⚠️ 처리된 QZA 파일이 없습니다.")
