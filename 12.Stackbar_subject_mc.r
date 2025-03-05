# StackbarExtended 패키지를 사용하여 metadata와 qiime exported biom tsv feature
# table 이용 Stacked Barplot 생성, PERMANOVA adonis2, Kruskal-Wallis rank sum
# test, FDR correction (BH method), 유의성 표기 추가

###############################################################################
# 0. 패키지 설치/로드 + 작업 디렉터리 설정
###############################################################################

install_and_load <- function(pkg, fromBioc = FALSE, fromGithub = FALSE, ghrepo = NULL, configure_vars = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing package:", pkg))
    if (!is.null(configure_vars)) {
      install.packages(pkg,
        repos = "http://cran.us.r-project.org",
        configure.vars = configure_vars
      )
    } else if (fromBioc) {
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
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# 필요한 패키지 설치+로드
install_and_load("phyloseq", fromBioc = TRUE) # phyloseq는 Bioconductor에서 설치
install_and_load("dplyr")
install_and_load("ggtext") # ggplot2 텍스트 레이블 서식 지정
install_and_load("DESeq2", fromBioc = TRUE) # DESeq2는 Bioconductor 패키지
install_and_load("ragg", configure_vars = "INCLUDE_DIR=/usr/include/freetype2 LIB_DIR=/usr/lib/x86_64-linux-gnu ZLIB_HOME=/usr")
install_and_load("StackbarExtended", fromGithub = TRUE, ghrepo = "ThibaultCuisiniere/StackbarExtended")
install_and_load("vegan") # PERMANOVA 분석을 위해 필요
install_and_load("compositions") # CLR 변환을 위해 필요
install_and_load("ggplot2") # ggsave() 등을 위해 필요

# 작업 폴더 설정
cwd <- "/Users/inseonghwang/OneDrive/Sparta_300"
# cwd <- "/mnt/d/OneDrive/Sparta_300"
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

max_len <- 7 # (Kingdom, Phylum, Class, Order, Family, Genus, Species) 최대 7개
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
# 6. PERMANOVA (Adonis2) + 개별 Taxa 검정 (Kruskal-Wallis) + 유의성 표기 + Genus 이탤릭 적용
###############################################################################

dist_matrix <- vegdist(t(otu_mat2), method = "euclidean") # 거리 행렬 생성

# PERMANOVA 실행
adonis2_result <- adonis2(
  dist_matrix ~ group,
  data = meta_df2,
  permutations = 999,
  by = "margin"
)

# PERMANOVA FDR 보정 적용
p_values <- adonis2_result$`Pr(>F)`
fdr_corrected_p <- p.adjust(p_values, method = "fdr")

# 개별 Taxa 차이를 평가하기 위한 Kruskal-Wallis 검정 수행
kruskal_p_values <- sapply(rownames(otu_mat2), function(taxa) {
  taxa_values <- otu_mat2[taxa, ] # 특정 박테리아의 모든 샘플 상대 빈도 가져오기
  kruskal.test(taxa_values ~ meta_df2$group)$p.value # Kruskal-Wallis 검정
})

# FDR 보정 적용
kruskal_fdr_corrected_p <- p.adjust(kruskal_p_values, method = "fdr")

# 유의성 기호 추가 (NA 값은 빈 문자열로 처리)
significance_levels <- ifelse(
  is.na(kruskal_fdr_corrected_p), "", # NA 값은 공백 처리
  ifelse(kruskal_fdr_corrected_p < 0.001, "***",
    ifelse(kruskal_fdr_corrected_p < 0.01, "**",
      ifelse(kruskal_fdr_corrected_p < 0.05, "*", "")
    )
  )
)

# 🔹 박테리아 Genus 이름을 이탤릭 적용 (Stacked Barplot 레전드 반영)
tax_df_clean$Genus <- paste0("*", tax_df_clean$Genus, "* ", significance_levels)

# phyloseq tax_table 업데이트
tax_mat_updated <- as.matrix(tax_df_clean)
TAX_updated <- phyloseq::tax_table(tax_mat_updated)

# phyloseq 객체 생성 (샘플별)
ps <- phyloseq(OTU2, TAX_updated, SAMP2)

# CSV 파일로 저장
adonis_out <- file.path(cwd, "07_collapse", "exported_barplot", "adonis2_results.csv")
dir.create(dirname(adonis_out), showWarnings = FALSE)
write.csv(as.data.frame(adonis2_result), adonis_out, row.names = TRUE)
message("PERMANOVA (Adonis2) 결과가 저장되었습니다: ", adonis_out)

###############################################################################
# 7. StackbarExtended → 레전드 크기 및 자동 여백 조정 후 Stacked Barplot 저장
###############################################################################

custom_theme <- ggplot2::theme(
  legend.text = ggtext::element_markdown(size = 12), # 🔹 이탤릭 적용된 텍스트 지원
  legend.title = ggplot2::element_text(size = 14, face = "bold"), # 🔹 레전드 제목 크기 확대
  legend.position = "right", # 🔹 레전드 위치 고정
  legend.box = "vertical", # 🔹 레전드를 세로 정렬하여 자동 확장
  plot.margin = ggplot2::margin(1, 2, 1, 1, unit = "cm") # 🔹 여백을 동적으로 조절
)

# Stacked Barplot 실행
res <- plot_microbiota(
  ps_object = ps,
  exp_group = "group",
  sample_name = "SampleID",
  main_level = "Phylum",
  sub_level = "Genus",
  differential_analysis = FALSE,
  sig_lab = TRUE, # 유의성 표시 추가
  hues = c("Blues", "Greens", "Oranges", "Purples", "Reds")
)

# 🔹 `res$plot`에 직접 테마 적용
res$plot <- res$plot + custom_theme

# 🔹 그래프 크기를 자동 조절하여 레전드가 완전히 포함되도록 설정
pdf_out <- file.path(cwd, "07_collapse", "exported_barplot", "StackbarExt_subject.pdf")

ggsave(
  filename = pdf_out, plot = res$plot,
  width = 14 + length(unique(tax_df_clean$Genus)) * 0.02, height = 10,
  limitsize = FALSE # 🔹 그래프 크기 제한 해제
)
message("Stacked Barplot이 저장되었습니다: ", pdf_out)
