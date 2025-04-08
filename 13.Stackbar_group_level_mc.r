# StackbarExtended 패키지를 사용하여 metadata와 qiime exported biom tsv feature
# table 이용 Stacked Barplot 생성, Friedman Test, FDR correction (BH method), 유의성 표기 추가

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
required_packages <- c("phyloseq", "dplyr", "ggtext", "vegan", "compositions", "ggplot2", "reshape2")
for (pkg in required_packages) {
    install_and_load(pkg)
}
install_and_load("StackbarExtended", fromGithub = TRUE, ghrepo = "ThibaultCuisiniere/StackbarExtended")

# 작업 폴더 설정
cwd <- "/Users/inseonghwang/OneDrive/Sparta_300"
setwd(cwd)

# 오류 발생 시 도움이 되는 디버깅 옵션 설정
options(error = function() {
    traceback(3)
    if (interactive()) recover()
})

###############################################################################
# 1. TSV (genus_20FNS-flt2-RF.tsv) 불러오기 → OTU 테이블 만들기
###############################################################################

tsv_path <- file.path(cwd, "07_collapse", "genus_20FNS-flt1-RF.tsv")
message(paste("파일 불러오기:", tsv_path))

# 파일 존재 여부 확인
if (!file.exists(tsv_path)) {
    stop("TSV 파일이 존재하지 않습니다. 경로를 확인하세요.")
}

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
message(paste("OTU 테이블 크기:", nrow(otu_mat), "OTUs x", ncol(otu_mat), "샘플"))

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

# 비어있는 taxa 채우기 (NA → "Unknown")
tax_df_clean[is.na(tax_df_clean) | tax_df_clean == ""] <- "Unknown"

# phyloseq tax_table 생성
tax_mat <- as.matrix(tax_df_clean)
TAX <- phyloseq::tax_table(tax_mat)

###############################################################################
# 3. 메타데이터 불러오기 → sample_data
###############################################################################

meta_file <- file.path(cwd, "metadata_20FNS.tsv")
# 파일 존재 여부 확인
if (!file.exists(meta_file)) {
    stop("메타데이터 파일이 존재하지 않습니다. 경로를 확인하세요.")
}

meta_df <- read.delim(
    meta_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    row.names = 1
)

# sample_name 인식 문제 해결을 위해 SampleID 열 추가
meta_df$SampleID <- rownames(meta_df)
message(paste("메타데이터 로드 완료:", nrow(meta_df), "샘플"))

###############################################################################
# 4. 공통 샘플만 필터링 (OTU vs 메타데이터)
###############################################################################

otu_samples <- colnames(otu_mat)
meta_samples <- rownames(meta_df)
common_samples <- intersect(otu_samples, meta_samples)
message(paste("OTU 테이블에 있는 샘플 수:", length(otu_samples)))
message(paste("메타데이터에 있는 샘플 수:", length(meta_samples)))
message(paste("공통 샘플 수:", length(common_samples)))

if (length(common_samples) == 0) {
    stop("OTU 테이블과 메타데이터 간 공통 샘플이 없습니다. 샘플 ID를 확인하세요.")
}

# 공통 샘플만 유지
otu_mat2 <- otu_mat[, common_samples, drop = FALSE]
meta_df2 <- meta_df[common_samples, , drop = FALSE]

# 0이 아닌 값을 가진 OTU만 유지 (선택적)
non_zero_rows <- rowSums(otu_mat2) > 0
otu_mat2 <- otu_mat2[non_zero_rows, , drop = FALSE]
message(paste("0이 아닌 값을 가진 OTU 수:", sum(non_zero_rows)))

# phyloseq 객체 생성을 위해 필요한 테이블 생성
OTU2 <- phyloseq::otu_table(otu_mat2, taxa_are_rows = TRUE)
SAMP2 <- sample_data(meta_df2)

###############################################################################
# 5. phyloseq 객체 생성
###############################################################################
ps <- phyloseq(OTU2, TAX[rownames(otu_mat2), ], SAMP2)
message("Phyloseq 객체 생성 완료.")

###############################################################################
# 6. 샘플을 DQ, MW, TC 그룹별로 평균 계산 후 phyloseq 객체 업데이트
# Friedman 테스트 및 FDR 보정으로 유의성 표시 추가
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
        message(paste(group, "그룹에서", length(sample_indices), "개 샘플 발견"))
    } else {
        warning(paste("그룹", group, "에 해당하는 샘플이 없습니다."))
    }
}

# phyloseq용 OTU 테이블 업데이트 (그룹 평균 적용)
OTU_grouped <- phyloseq::otu_table(as.matrix(otu_grouped), taxa_are_rows = TRUE)

# 새로운 그룹별 메타데이터 생성
meta_grouped <- data.frame(Group = group_labels, row.names = group_labels)
SAMP_grouped <- phyloseq::sample_data(meta_grouped)

# Friedman 검정 및 FDR 보정 (BH 방법) 수행
significance_levels <- rep("", nrow(tax_df_clean[rownames(otu_mat2), ]))

# 메타데이터에서 피험자 ID 추출
# 샘플 이름 형식이 'Group_SubjectID'로 되어 있다고 가정
sample_names <- colnames(otu_mat2)
subject_ids <- sub("^[^_]+_", "", sample_names) # 'Group_' 부분 제거하여 SubjectID만 추출
subject_groups <- sub("_.*$", "", sample_names) # 그룹 정보 추출 (DQ, MW, TC)

# 피험자별 그룹 수 확인
subject_group_counts <- table(subject_ids)
complete_subjects <- names(subject_group_counts[subject_group_counts == length(group_labels)])
message(paste("모든 그룹에 데이터가 있는 피험자 수:", length(complete_subjects)))

if (length(complete_subjects) < 3) {
    warning("Friedman 검정을 위한 충분한 데이터가 없습니다. 최소 3명의 피험자가 모든 그룹에 데이터를 가져야 합니다.")
}

# 각 taxa에 대해 Friedman 검정 수행
p_values <- numeric(nrow(otu_mat2))
names(p_values) <- rownames(otu_mat2)

for (i in 1:nrow(otu_mat2)) {
    # 데이터프레임 생성: 피험자 ID, 그룹, 해당 taxa의 abundance
    test_data <- data.frame(
        subject = subject_ids,
        group = subject_groups,
        abundance = as.numeric(otu_mat2[i, ])
    )

    # 각 피험자가 모든 그룹(DQ, MW, TC)의 데이터가 있는지 확인
    valid_subjects <- intersect(unique(test_data$subject), complete_subjects)

    if (length(valid_subjects) >= 3) { # 최소 3명 이상의 완전한 데이터가 있는 피험자 필요
        # 완전한 데이터가 있는 피험자만 선택
        complete_data <- test_data[test_data$subject %in% valid_subjects, ]

        # 데이터를 넓은 형식(wide format)으로 변환
        wide_data <- reshape2::dcast(complete_data, subject ~ group, value.var = "abundance")

        # 첫 번째 열(subject)을 제외한 그룹 데이터만 추출
        group_data <- as.matrix(wide_data[, -1])

        # Friedman 검정 수행
        tryCatch(
            {
                friedman_result <- friedman.test(group_data)
                p_values[i] <- friedman_result$p.value
            },
            error = function(e) {
                warning(paste("Friedman 검정 오류 (row", i, "):", e$message))
                p_values[i] <- 1.0
            }
        )
    } else {
        p_values[i] <- 1.0 # 충분한 데이터가 없는 경우 유의하지 않음으로 설정
    }
}

# BH 방법으로 p-value 보정
adjusted_p <- p.adjust(p_values, method = "BH")

# 보정된 p-value에 따라 별표 재할당
significance_levels <- rep("", length(adjusted_p))
significance_levels[adjusted_p < 0.05] <- "*"
significance_levels[adjusted_p < 0.01] <- "**"
significance_levels[adjusted_p < 0.001] <- "***"

# 유의한 차이를 보이는 taxa 수 출력
message(paste("Friedman 검정 후 유의한 차이를 보이는 taxa (p < 0.05):", sum(adjusted_p < 0.05)))
message(paste("유의한 차이를 보이는 taxa (p < 0.01):", sum(adjusted_p < 0.01)))
message(paste("유의한 차이를 보이는 taxa (p < 0.001):", sum(adjusted_p < 0.001)))

# Genus 이름에 유의성 표시 추가 후 TAX_updated 생성
tax_df_formatted <- tax_df_clean[rownames(otu_mat2), ]
tax_df_formatted$Genus <- paste0("*", tax_df_formatted$Genus, "*", significance_levels)

# tax_table 업데이트
tax_mat_updated <- as.matrix(tax_df_formatted)
TAX_updated <- phyloseq::tax_table(tax_mat_updated)

# 새로운 phyloseq 객체 생성 (그룹 평균 적용)
ps_grouped <- phyloseq(OTU_grouped, TAX_updated, SAMP_grouped)

# `ps_grouped` 객체 생성 확인
if (!exists("ps_grouped")) {
    stop("Error: ps_grouped object was not created correctly.")
}
message("그룹별 평균 기반 Phyloseq 객체 생성 완료 (Genus 이탤릭 적용됨, Friedman test 유의성 표시 추가).")

###############################################################################
# 7. StackbarExtended → 레전드 크기 및 자동 여백 조정 후 Stacked Barplot 저장
###############################################################################

# 출력 디렉터리 생성
output_dir <- file.path(cwd, "07_collapse", "exported_barplot")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

custom_theme <- ggplot2::theme(
    legend.text = ggtext::element_markdown(size = 12), # 이탤릭 적용된 텍스트 지원
    legend.title = ggplot2::element_text(size = 14, face = "bold"), # 레전드 제목 크기 확대
    legend.position = "right", # 레전드 위치 고정
    legend.box = "vertical", # 레전드를 세로로 정렬하여 자동 확장
    plot.margin = ggplot2::margin(1, 2, 1, 1, unit = "cm"), # 여백을 동적으로 조절
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 12), # X축 텍스트 회전 및 크기 조정
    axis.title = ggplot2::element_text(size = 14, face = "bold") # 축 제목 크기 및 스타일 조정
)

# Stacked Barplot 실행
tryCatch(
    {
        res <- plot_microbiota(
            ps_object = ps_grouped, # 그룹별 평균이 반영된 phyloseq 객체 사용
            exp_group = "Group", # 그룹 정보 반영
            sample_name = "Group", # 그룹 이름 사용
            main_level = "Phylum",
            sub_level = "Genus",
            n_phy = 7, # 상위 7개 Phylum만 표시
            differential_analysis = FALSE,
            sig_lab = TRUE, # 유의성 표시 추가
            #   hues = c("Blues", "Purples", "Greens", "Reds", "Oranges", "Dark2", "Paired"), # 색상 팔레트 지정
            hues = c("Blues", "Purples", "Paired", "Reds", "Greens", "Oranges", "Dark2"), # 색상 팔레트 지정
        )

        # `res$plot`에 직접 테마 적용
        res$plot <- res$plot + custom_theme

        # 그래프 크기를 자동 조절하여 레전드가 완전히 포함되도록 설정
        pdf_out <- file.path(output_dir, "StackbarExt_group.pdf")
        png_out <- file.path(output_dir, "StackbarExt_group.png")

        # PDF 저장
        ggsave(
            filename = pdf_out, plot = res$plot,
            width = 2 + length(unique(tax_df_clean$Genus)) * 0.02, height = 10,
            limitsize = FALSE # 그래프 크기 제한 해제
        )

        # PNG도 함께 저장 (웹 보기 용이)
        ggsave(
            filename = png_out, plot = res$plot,
            width = 12, height = 10, dpi = 300,
            limitsize = FALSE
        )

        message("그룹별 평균 기반 Stacked Barplot이 저장되었습니다:")
        message(paste(" - PDF:", pdf_out))
        message(paste(" - PNG:", png_out))
    },
    error = function(e) {
        message("Stacked Barplot 생성 중 오류 발생:", e$message)
    }
)

###############################################################################
# 8. 각 그룹별 genus 레벨의 퍼센트 테이블 저장
###############################################################################

# genus 레벨 데이터 추출
save_genus_percentages_by_group <- function(ps_obj, output_dir) {
    # 디렉터리 생성 (없는 경우)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # Genus 레벨 데이터 추출
    genus_data <- tax_glom(ps_obj, taxrank = "Genus")

    # OTU 테이블 추출
    otu_table_df <- as.data.frame(otu_table(genus_data))

    # Taxonomy 테이블 추출
    tax_table_df <- as.data.frame(tax_table(genus_data))

    # Genus 이름 정리 (별표 등 제거)
    clean_genus_names <- gsub("\\*|\\s+\\*+$", "", tax_table_df$Genus)

    # 결과 데이터프레임 생성
    result_df <- data.frame(
        Genus = clean_genus_names,
        Phylum = tax_table_df$Phylum
    )

    # 각 그룹별 퍼센트 추가
    for (group in colnames(otu_table_df)) {
        # 퍼센트로 변환 (100 곱하기)
        result_df[[group]] <- otu_table_df[, group] * 100
    }

    # 각 그룹별로 정렬된 파일 생성
    for (group in colnames(otu_table_df)) {
        # 현재 그룹에 대한 데이터 정렬
        group_df <- result_df[order(result_df[[group]], decreasing = TRUE), ]

        # CSV 파일로 저장
        output_file <- file.path(output_dir, paste0("genus_percentage_", group, ".csv"))
        write.csv(group_df, file = output_file, row.names = FALSE)
        message(paste0(group, " 그룹의 genus 퍼센트 테이블이 저장되었습니다: ", output_file))
    }

    # 전체 데이터 저장 (모든 그룹 포함)
    all_groups_file <- file.path(output_dir, "genus_percentage_all_groups.csv")
    write.csv(result_df, file = all_groups_file, row.names = FALSE)
    message("모든 그룹의 genus 퍼센트 테이블이 저장되었습니다: ", all_groups_file)

    return(result_df)
}

# 함수 실행
output_dir <- file.path(cwd, "07_collapse", "exported_barplot", "genus_tables")
genus_percentages <- save_genus_percentages_by_group(ps_grouped, output_dir)

# 추가로 Phylum별로 그룹화된 테이블 생성
phylum_grouped_file <- file.path(output_dir, "genus_percentage_by_phylum.csv")
genus_by_phylum <- genus_percentages %>%
    arrange(Phylum, desc(DQ)) %>%
    select(Genus, Phylum, everything())
write.csv(genus_by_phylum, file = phylum_grouped_file, row.names = FALSE)
message("Phylum별로 정렬된 genus 퍼센트 테이블이 저장되었습니다: ", phylum_grouped_file)

# 스크립트 종료 시간 기록
end_time <- Sys.time()
total_time <- difftime(end_time, Sys.time(), units = "mins")
message(paste("분석 완료! 총 실행 시간:", round(total_time, 2), "분"))
