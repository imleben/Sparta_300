# StackbarExtended 패키지를 사용하여 metadata와 qiime exported biom tsv feature
# table 이용 Stacked Barplot 생성, PERMANOVA adonis2, Kruskal-Wallis rank sum
# test, FDR correction (BH method), 유의성 표기 추가

###############################################################################
# 0. 패키지 설치/로드 + 작업 디렉터리 설정
###############################################################################

install_and_load <- function(pkg, fromBioc = FALSE, fromGithub = FALSE, ghrepo = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (fromBioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "http://cran.us.r-project.org")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else if (fromGithub) {
      if (!requireNamespace("devtools", quietly = TRUE)) {
        install.packages("devtools", repos = "http://cran.us.r-project.org")
      }
      devtools::install_github(ghrepo)
    } else {
      install.packages(pkg, repos = "http://cran.us.r-project.org")
    }
  }
  library(pkg, character.only = TRUE)
}

# 필요한 패키지 설치+로드
install_and_load("phyloseq")
install_and_load("dplyr")
install_and_load("ggtext") # ggplot2 텍스트 레이블 서식 지정
install_and_load("StackbarExtended", fromGithub = TRUE, ghrepo = "ThibaultCuisiniere/StackbarExtended")
install_and_load("vegan") # PERMANOVA 분석을 위해 필요
install_and_load("compositions") # CLR 변환을 위해 필요
install_and_load("ggplot2") # ggsave() 위해 필요

# 작업 폴더 설정
cwd <- "/Users/inseonghwang/OneDrive/Sparta_300"
setwd(cwd)

###############################################################################
# 1. TSV (genus_20FNS-flt2-RF.tsv) 불러오기 → OTU 테이블 만들기
###############################################################################

tsv_path <- file.path(cwd, "07_collapse", "genus_20FNS-flt1-RF.tsv")
df <- read.delim(
  tsv_path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 첫 열 이름 재설정 (예: '#OTU ID' → 'taxonomy')
colnames(df)[1] <- "taxonomy"
rownames(df) <- df$taxonomy

# 첫 열 제외 → OTU(행) vs 샘플(열) 행렬
otu_mat <- as.matrix(df[, -1])
mode(otu_mat) <- "numeric" # 상대 빈도(또는 비율)로 가정

# phyloseq용 OTU 테이블 생성
OTU <- phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)

###############################################################################
# 2. Taxonomy 문자열 파싱 → tax_table 생성
###############################################################################

tax_strings <- rownames(df)
split_tax <- strsplit(tax_strings, ";\\s*")

max_len <- 7 # (Domain, Phylum, Class, Order, Family, Genus, Species) 최대 7개
split_tax_pad <- lapply(split_tax, function(x) {
  length(x) <- max_len
  x
})

tax_df <- data.frame(
  do.call(rbind, split_tax_pad),
  stringsAsFactors = FALSE
)
colnames(tax_df) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(tax_df) <- rownames(df)

# d__, p__ 접두사 제거
remove_prefix <- function(x) sub("^[a-z]__", "", x)
tax_df_clean <- as.data.frame(apply(tax_df, 2, remove_prefix), stringsAsFactors = FALSE)
colnames(tax_df_clean) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# phyloseq tax_table 생성
tax_mat <- as.matrix(tax_df_clean)
TAX <- phyloseq::tax_table(tax_mat)

###############################################################################
# 3. 메타데이터 불러오기 → sample_data
###############################################################################

meta_df <- read.delim(
  "metadata_20FNS.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  row.names = 1
)

# sample_name 인식 문제 해결을 위해 SampleID 열 추가
meta_df$SampleID <- rownames(meta_df)

###############################################################################
# 4. 공통 샘플만 필터링 (OTU vs 메타데이터)
###############################################################################

otu_samples <- colnames(otu_mat)
meta_samples <- rownames(meta_df)
common_samples <- intersect(otu_samples, meta_samples)

otu_mat2 <- otu_mat[, common_samples, drop = FALSE]
meta_df2 <- meta_df[common_samples, , drop = FALSE]

# phyloseq 객체 생성을 위해 필요한 테이블 생성
OTU2 <- phyloseq::otu_table(otu_mat2, taxa_are_rows = TRUE)
SAMP2 <- sample_data(meta_df2)

###############################################################################
# 5. phyloseq 객체 생성
###############################################################################
ps <- phyloseq(OTU2, TAX, SAMP2)
message("Phyloseq 객체 생성 완료.")

###############################################################################
# 6. 샘플을 DQ, MW, TC 그룹별로 평균 계산 후 phyloseq 객체 업데이트 (레전드 이탤릭 적용)
###############################################################################

# (1) DQ, MW, TC 그룹별로 샘플 분류
group_labels <- c("DQ", "MW", "TC")

# 그룹별 평균을 계산할 데이터프레임 생성
otu_grouped <- data.frame(matrix(ncol = length(group_labels), nrow = nrow(otu_mat2)))
colnames(otu_grouped) <- group_labels
rownames(otu_grouped) <- rownames(otu_mat2)

# 각 그룹에 해당하는 샘플을 찾아 평균 계산
for (group in group_labels) {
  sample_indices <- grep(paste0("^", group, "_"), colnames(otu_mat2)) # DQ_, MW_, TC_로 시작하는 샘플 찾기
  if (length(sample_indices) > 0) {
    otu_grouped[, group] <- rowMeans(otu_mat2[, sample_indices, drop = FALSE]) # 평균 계산
  } else {
    message("Warning: No samples found for group ", group)
  }
}

# 🔹 phyloseq용 OTU 테이블 업데이트 (그룹 평균 적용)
OTU_grouped <- phyloseq::otu_table(as.matrix(otu_grouped), taxa_are_rows = TRUE)

# 새로운 그룹별 메타데이터 생성
meta_grouped <- data.frame(Group = group_labels, row.names = group_labels)
SAMP_grouped <- phyloseq::sample_data(meta_grouped)

# 🔹 phyloseq tax_table 업데이트 (Genus 이름을 이탤릭 적용)
if (!exists("TAX_updated")) {
  stop("Error: TAX_updated object was not created correctly.")
}

# Genus 이름을 이탤릭으로 변경
tax_df_clean$Genus <- paste0("*", tax_df_clean$Genus, "* ", significance_levels)

# tax_table 업데이트
tax_mat_updated <- as.matrix(tax_df_clean)
TAX_updated <- phyloseq::tax_table(tax_mat_updated)

# 🔹 새로운 phyloseq 객체 생성 (그룹 평균 적용)
ps_grouped <- phyloseq(OTU_grouped, TAX_updated, SAMP_grouped)

# 🔹 `ps_grouped` 객체 생성 확인
if (!exists("ps_grouped")) {
  stop("Error: ps_grouped object was not created correctly.")
}
message("그룹별 평균 기반 Phyloseq 객체 생성 완료 (Genus 이탤릭 적용됨).")

###############################################################################
# 7. StackbarExtended → 레전드 크기 및 자동 여백 조정 후 Stacked Barplot 저장
###############################################################################

custom_theme <- ggplot2::theme(
  legend.text = ggtext::element_markdown(size = 12), # 🔹 이탤릭 적용된 텍스트 지원
  legend.title = ggplot2::element_text(size = 14, face = "bold"), # 🔹 레전드 제목 크기 확대
  legend.position = "right", # 🔹 레전드 위치 고정
  legend.box = "vertical", # 🔹 레전드를 세로로 정렬하여 자동 확장
  plot.margin = ggplot2::margin(1, 2, 1, 1, unit = "cm") # 🔹 여백을 동적으로 조절
)

# Stacked Barplot 실행
res <- plot_microbiota(
  ps_object = ps_grouped, # 그룹별 평균이 반영된 phyloseq 객체 사용
  exp_group = "Group", # 그룹 정보 반영
  sample_name = "Group", # 그룹 이름 사용
  main_level = "Phylum",
  sub_level = "Genus",
  differential_analysis = FALSE,
  sig_lab = TRUE, # 유의성 표시 포함
  hues = c("Blues", "Greens", "Oranges", "Purples", "Reds")
)

# 🔹 `res$plot`에 직접 테마 적용
res$plot <- res$plot + custom_theme

# 🔹 그래프 크기를 자동 조절하여 레전드가 완전히 포함되도록 설정
pdf_out <- file.path(cwd, "07_collapse", "exported_barplot", "StackbarExt_group.pdf")
ggsave(
  filename = pdf_out, plot = res$plot,
  width = 2 + length(unique(tax_df_clean$Genus)) * 0.02, height = 10,
  limitsize = FALSE # 🔹 그래프 크기 제한 해제
)
message("그룹별 평균 기반 Stacked Barplot이 저장되었습니다: ", pdf_out)
