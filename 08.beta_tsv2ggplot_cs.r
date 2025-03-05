# PCoA meta combined complete set tsv 파일을 읽어 PCoA 결과를 시각화하고 PERMANOVA를 수행하는 R 스크립트
# (추가: PERMANOVA adonis2 결과 중 DQ vs MW, MW vs TC, TC vs DQ 그룹별 비교 p‑value 저장)

# 1. 패키지 로드
required_packages <- c("ggplot2", "dplyr", "ellipse", "vegan")
for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, repos = "https://cran.rstudio.com/")
        library(pkg, character.only = TRUE)
    }
}
print("All required packages are successfully installed and loaded.")

# 2. PCoA 결과 파일이 들어있는 폴더 (이미 메타데이터가 포함된 TSV)
pcoa_dir <- "/Users/inseonghwang/OneDrive/Sparta_300/05_diversity_20FNS/beta_meta_combined/pcoa2tsv_all"
pcoa_files <- list.files(path = pcoa_dir, pattern = "\\.tsv$", full.names = TRUE)

# 출력 폴더 경로 (여기서는 pcoa_dir와 동일한 폴더 사용)
output_dir <- pcoa_dir
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# 3. PERMANOVA 결과를 저장할 객체
permanova_results <- data.frame(
    File = character(),
    Overall_P_Value = character(),
    P_DQ_vs_MW = character(),
    P_MW_vs_TC = character(),
    P_TC_vs_DQ = character(),
    stringsAsFactors = FALSE
)

# 4. 각 PCoA 파일 순회
for (pcoa_file in pcoa_files) {
    # (a) PCoA TSV 읽기 (이미 group, subject, PC1~3 등이 포함)
    pcoa_data <- read.table(pcoa_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    # group 컬럼이 있는지 확인
    if (!"group" %in% colnames(pcoa_data)) {
        stop(sprintf("Error: 'group' column not found in %s", basename(pcoa_file)))
    }
    # subject 컬럼이 있는지 확인 (PERMANOVA strata 용)
    if (!"subject" %in% colnames(pcoa_data)) {
        stop(sprintf("Error: 'subject' column not found in %s", basename(pcoa_file)))
    }

    # -------------------
    # (1) ggplot으로 2D 시각화 (PC1 vs PC2)
    # -------------------
    p <- ggplot(pcoa_data, aes(x = PC1, y = PC2, color = group, fill = group)) +
        geom_point(shape = 21, color = "black", size = 5, alpha = 0.7) +
        stat_ellipse(geom = "polygon", alpha = 0.2, linetype = "blank", level = 0.95) +
        scale_fill_manual(values = c("#ec2424", "#08b008", "#1d1de9")) +
        scale_x_continuous(labels = function(x) sprintf("%.1f", x)) +
        scale_y_continuous(labels = function(x) sprintf("%.1f", x)) +
        theme_minimal() +
        theme(
            axis.title.x = element_text(size = 28, margin = margin(t = 15)),
            axis.title.y = element_text(size = 28, margin = margin(r = 5)),
            axis.text = element_text(size = 24),
            plot.title = element_text(size = 14, margin = margin(b = 10)),
            legend.position = "right",
            legend.text = element_text(size = 26),
            legend.title = element_blank(),
            legend.key.height = unit(1.5, "cm"),
            legend.spacing.y = unit(0.5, "cm"),
            # 아래 x축과 왼쪽 y축을 검정색 실선으로 설정:
            axis.line.x = element_line(color = "black"),
            axis.line.y = element_line(color = "black")
        ) +
        labs(
            x = "PC1",
            y = "PC2",
            title = gsub("_pcoa_results.tsv", "", basename(pcoa_file))
        ) +
        guides(
            color = guide_legend(override.aes = list(size = 6)),
            fill  = guide_legend(override.aes = list(size = 6))
        )

    # 그래프 저장 (PDF)
    pdf_file <- file.path(output_dir, gsub(".tsv$", ".pdf", basename(pcoa_file)))
    ggsave(filename = pdf_file, width = 10, height = 10, plot = p)

    # -------------------
    # (2) PERMANOVA (PC1, PC2, PC3) - 전체 그룹 비교
    # -------------------
    needed_cols <- c("PC1", "PC2", "PC3")
    missing_cols <- setdiff(needed_cols, colnames(pcoa_data))
    if (length(missing_cols) > 0) {
        warning(sprintf(
            "[%s] Missing columns for PERMANOVA: %s -> Skipping",
            basename(pcoa_file),
            paste(missing_cols, collapse = ", ")
        ))
        next
    }
    if (nrow(pcoa_data) < 3) {
        warning(sprintf("[%s] Too few samples (<3) for PERMANOVA -> Skipping", basename(pcoa_file)))
        next
    }

    # 거리 행렬 (유클리드)
    dist_matrix <- dist(pcoa_data[, needed_cols], method = "euclidean")
    adonis_result <- adonis2(
        dist_matrix ~ group,
        data = pcoa_data,
        permutations = 9999,
        strata = pcoa_data$subject
    )
    raw_overall_p <- adonis_result$`Pr(>F)`[1]
    overall_p_str <- ifelse(raw_overall_p < 0.001, "<0.001", formatC(raw_overall_p, format = "f", digits = 4))

    # -------------------
    # (3) Pairwise PERMANOVA comparisons: DQ vs MW, MW vs TC, TC vs DQ
    # -------------------
    comparisons <- list(
        "DQ_vs_MW" = c("DQ", "MW"),
        "MW_vs_TC" = c("MW", "TC"),
        "TC_vs_DQ" = c("TC", "DQ")
    )
    pairwise_pvals <- list()
    for (comp_name in names(comparisons)) {
        comp_groups <- comparisons[[comp_name]]
        subset_data <- pcoa_data[pcoa_data$group %in% comp_groups, ]
        if (length(unique(subset_data$group)) < 2) {
            pairwise_pvals[[comp_name]] <- NA
            next
        }
        dist_subset <- dist(subset_data[, needed_cols], method = "euclidean")
        adonis_pair <- adonis2(
            dist_subset ~ group,
            data = subset_data,
            permutations = 9999
        )
        raw_pair_p <- adonis_pair$`Pr(>F)`[1]
        pairwise_pvals[[comp_name]] <- ifelse(raw_pair_p < 0.001, "<0.001", formatC(raw_pair_p, format = "f", digits = 4))
    }

    # -------------------
    # (4) PERMANOVA 결과 저장 (전체 및 pairwise)
    # -------------------
    permanova_results <- rbind(
        permanova_results,
        data.frame(
            File = basename(pcoa_file),
            Overall_P_Value = overall_p_str,
            P_DQ_vs_MW = pairwise_pvals[["DQ_vs_MW"]],
            P_MW_vs_TC = pairwise_pvals[["MW_vs_TC"]],
            P_TC_vs_DQ = pairwise_pvals[["TC_vs_DQ"]],
            stringsAsFactors = FALSE
        )
    )
}

# -------------------
# (5) 모든 PCoA 파일 처리 후, PERMANOVA 결과 CSV로 저장
# -------------------
csv_file <- file.path(output_dir, "permanova_results.csv")
write.csv(permanova_results, csv_file, row.names = FALSE)
message(sprintf("PERMANOVA results saved to: %s", csv_file))
