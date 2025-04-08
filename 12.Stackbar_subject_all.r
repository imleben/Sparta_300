# StackbarExtended 패키지를 사용하여 metadata와 qiime exported biom tsv feature
# table 이용 Stacked Barplot 생성

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
install_and_load("ggplot2") # ggsave() 등을 위해 필요

# 작업 폴더 설정
cwd <- "/Users/inseonghwang/OneDrive/Sparta_300"
# cwd <- "/mnt/d/OneDrive/Sparta_300"
setwd(cwd)

###############################################################################
# 1. TSV (genus_20FNS-flt1-RF.tsv) 불러오기 → OTU 테이블 만들기
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

# 비어있는 분류군 값을 "Unknown"으로 대체
tax_df_clean[tax_df_clean == ""] <- "Unknown"

# Genus 이름을 이탤릭으로 적용
tax_df_clean$Genus <- paste0("*", tax_df_clean$Genus, "*")

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
# 6. StackbarExtended → 레전드 크기 및 자동 여백 조정 후 Stacked Barplot 저장
###############################################################################

custom_theme <- ggplot2::theme(
    legend.text = ggtext::element_markdown(size = 12), # 이탤릭 적용된 텍스트 지원
    legend.title = ggplot2::element_text(size = 14, face = "bold"), # 레전드 제목 크기 확대
    legend.position = "right", # 레전드 위치 고정
    legend.box = "vertical", # 레전드를 세로 정렬하여 자동 확장
    plot.margin = ggplot2::margin(1, 2, 1, 1, unit = "cm") # 여백을 동적으로 조절
)

# Stacked Barplot 실행
res <- plot_microbiota(
    ps_object = ps,
    exp_group = "group",
    sample_name = "SampleID",
    main_level = "Phylum",
    sub_level = "Genus",
    n_phy = 7, # 상위 7개 Phylum만 표시
    differential_analysis = FALSE,
    sig_lab = FALSE, # 유의성 표시 비활성화
    hues = c("Blues", "Purples", "Paired", "Reds", "Greens", "Oranges", "Dark2"), # 색상 팔레트 지정
)

# `res$plot`에 직접 테마 적용
res$plot <- res$plot + custom_theme

# 그래프 크기를 자동 조절하여 레전드가 완전히 포함되도록 설정
pdf_out <- file.path(cwd, "07_collapse", "exported_barplot", "StackbarExt_subject_all.pdf")
dir.create(dirname(pdf_out), showWarnings = FALSE, recursive = TRUE)

ggsave(
    filename = pdf_out, plot = res$plot,
    width = 14 + length(unique(tax_df_clean$Genus)) * 0.02, height = 10,
    limitsize = FALSE # 그래프 크기 제한 해제
)
message("Stacked Barplot이 저장되었습니다: ", pdf_out)
