# PCOA QZA 파일을 TSV로 변환하고 메타데이터와 병합하는 과정을 반복하는 스크립트.

# 1. 패키지 로드
required_packages <- c("qiime2R", "dplyr", "readr")

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("jbisanz/qiime2R")

# 설치되지 않은 패키지 확인 및 설치
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cran.rstudio.com/")
    library(pkg, character.only = TRUE)
  }
}

# 모든 패키지가 로드되었는지 확인
print("All required packages are successfully installed and loaded.")

# 2. 파일 경로 설정
qza_folder <- "/users/inseonghwang/onedrive/Sparta_300/05_diversity_20FNS"
metadata_path <- "/users/inseonghwang/onedrive/Sparta_300/metadata_20FNS.tsv"

tsv_output_folder <- "/users/inseonghwang/onedrive/Sparta_300/05_diversity_20FNS/beta_meta_combined"
merged_output_folder <- "/users/inseonghwang/onedrive/Sparta_300/05_diversity_20FNS/beta_meta_combined/pcoa2tsv_all"

# 3. 폴더 생성 (없으면)
if (!dir.exists(tsv_output_folder)) {
  dir.create(tsv_output_folder, recursive = TRUE)
}
if (!dir.exists(merged_output_folder)) {
  dir.create(merged_output_folder, recursive = TRUE)
}

# 4. 메타데이터 읽기
#    첫 번째 컬럼이 sample.id라고 가정 (다른 이름이면 강제로 바꾸기)
metadata <- read.delim(metadata_path, sep = "\t", stringsAsFactors = FALSE)
colnames(metadata)[1] <- "sample.id" # 필요하면 강제로 변경

# 5. QZA 파일 목록
qza_files <- list.files(
  path = qza_folder,
  pattern = "*_pcoa_results.qza",
  full.names = TRUE
)

# 6. 변환 & 병합 처리
for (qza_file in qza_files) {
  # (a) QZA 읽기
  pcoa_qza <- read_qza(qza_file)

  # (b) 출력될 TSV 파일 이름 (.qza -> .tsv)
  tsv_filename <- sub("\\.qza$", ".tsv", basename(qza_file))
  tsv_path <- file.path(tsv_output_folder, tsv_filename)

  # (c) pcoa_qza$data 구조 확인 후 TSV로 저장
  if (is.data.frame(pcoa_qza$data)) {
    # 전체가 df인 경우
    pcoa_df <- pcoa_qza$data
    write_tsv(pcoa_df, tsv_path)
  } else if (is.list(pcoa_qza$data) && "Vectors" %in% names(pcoa_qza$data)) {
    pcoa_df <- pcoa_qza$data$Vectors
    write_tsv(pcoa_df, tsv_path)
  } else {
    stop("pcoa_qza$data is not a data frame or does not contain 'Vectors'")
  }

  # (d) 병합 전 컬럼명 맞추기
  #     (예: PCoA 결과에는 "SampleID"라 되어 있을 수 있으므로 변경)
  colnames(pcoa_df)[colnames(pcoa_df) == "SampleID"] <- "sample.id"

  # (e) 메타데이터 병합 (inner join)
  merged_df <- pcoa_df %>%
    inner_join(metadata, by = "sample.id")
  # inner_join 이므로 PCoA TSV에 없는 sample.id는 제외

  # (f) 최종 병합 TSV 파일 이름 (.tsv -> _meta.tsv)
  merged_filename <- sub("\\.tsv$", "_meta_all.tsv", tsv_filename)
  merged_path <- file.path(merged_output_folder, merged_filename)

  # (g) 병합 결과 저장
  write_tsv(merged_df, merged_path)

  # (h) 진행 상황 출력
  message(sprintf(
    "[Done] %s -> %s & %s (Rows in merged: %d)",
    basename(qza_file),
    tsv_filename,
    merged_filename,
    nrow(merged_df)
  ))
}
