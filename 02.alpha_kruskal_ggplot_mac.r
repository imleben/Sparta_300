# Alpha diversity 지표에 대해 Kruskal-Wallis 검정을 수행하고, 그래프를 생성하는 스크립트

# 필요한 패키지 목록
required_packages <- c("ggplot2", "dplyr", "tidyr", "ggsignif", "tibble", "scales")

# 설치되지 않은 패키지 확인 및 설치
for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, repos = "https://cran.rstudio.com/")
        library(pkg, character.only = TRUE)
    }
}

# 모든 패키지가 로드되었는지 확인
print("All required packages are successfully installed and loaded.")

# 데이터 읽기
file_path <- "/Users/inseonghwang/onedrive/Sparta_300/05_diversity_20FNS/alpha_meta_combined/alpha_meta_combined.tsv"
data <- read.csv(file_path, sep = "\t")
data$group <- as.factor(data$group)

# 출력 경로
output_dir <- "/Users/inseonghwang/OneDrive/Sparta_300/05_diversity_20FNS/alpha_meta_combined/alpha_kruskal"

# 그래프를 생성하고 저장하는 함수 정의
plot_and_save <- function(metric, data, output_dir) {
    # Kruskal-Wallis 검정
    kruskal_test <- kruskal.test(reformulate("group", metric), data = data)
    kruskal_p <- kruskal_test$p.value

    # Wilcoxon pairwise 검정
    pairwise_p <- pairwise.wilcox.test(
        data[[metric]],
        data$group,
        p.adjust.method = "BH",
        exact = FALSE
    )

    # pairwise 결과 처리
    pairwise_results <- as.data.frame(as.table(pairwise_p$p.value)) %>%
        filter(!is.na(Freq)) %>%
        rename(Group1 = Var1, Group2 = Var2, p.value = Freq) %>%
        mutate(signif = case_when(
            p.value < 0.001 ~ "***",
            p.value < 0.01 ~ "**",
            p.value < 0.05 ~ "*",
            TRUE ~ ""
        )) %>%
        filter(signif != "")

    # y축 최대값, 최소값, 범위 계산
    max_value <- max(data[[metric]], na.rm = TRUE)
    min_value <- min(data[[metric]], na.rm = TRUE)
    range_value <- max_value - min_value

    # y축 여백 설정
    lower_margin <- min_value - range_value * 0.3 # 아래쪽 여백 (최대 30%)
    upper_margin <- max_value + range_value * 0.5 # 위쪽 여백 (최대 50%)

    # 브라켓 위치를 y축 범위에 상대적으로 설정
    y_bracket_positions <- max_value + seq(0.05, 0.15, length.out = 3) * range_value # 상대적 위치

    # p-value 텍스트 위치 계산
    p_value_position <- min_value - range_value * 0.2 # min_value 기준 상대 위치

    # p-value를 지수형으로 변환하여 위첨자로 표시하고 소수점 자릿수를 2자리로 제한
    format_p_value <- function(p) {
        exponent <- floor(log10(p))
        base <- format(p / 10^exponent, digits = 2, nsmall = 2)
        return(bquote(italic(p) == .(base) %*% 10^.(exponent)))
    }

    # 그래프 생성
    p <- ggplot(data, aes(x = group, y = .data[[metric]], fill = group)) +
        geom_violin(alpha = 0.5, color = NA, scale = "width", trim = FALSE) +
        geom_boxplot(width = 0.4, color = "black", alpha = 0.8, outlier.shape = NA) +
        geom_jitter(width = 0.2, size = 4, alpha = 0.3) +
        geom_signif(
            comparisons = list(c("DQ", "MW"), c("MW", "TC"), c("DQ", "TC")),
            map_signif_level = TRUE,
            y_position = y_bracket_positions, # 상대적 브라켓 위치
            step_increase = 0.1,
            tip_length = 0.0,
            textsize = 10
        ) +
        scale_y_continuous(
            limits = c(lower_margin, upper_margin), # y축 최소/최대값 설정
            expand = expansion(mult = c(0, 0.05)),
            breaks = scales::pretty_breaks(n = 5)
        ) +
        theme_minimal(base_size = 20) +
        theme(
            axis.title = element_text(size = 30),
            axis.text = element_text(size = 20, color = "black"),
            axis.line = element_line(color = "black"),
            panel.grid.major.y = element_line(color = "white"),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            legend.position = "none"
        ) +
        labs(
            title = metric,
            x = "",
            y = "Index"
        ) +
        annotate(
            "text",
            x = 1.5,
            y = p_value_position, # p-value 위치를 상대적으로 계산
            label = format_p_value(kruskal_p),
            size = 8,
            hjust = 0.2,
            vjust = 0.8
        )

    # 그래프 출력
    print(p)

    # PDF로 저장
    output_file <- file.path(output_dir, paste0(metric, "_plot.pdf"))
    ggsave(
        filename = output_file,
        plot = p,
        width = 5,
        height = 10,
        device = "pdf"
    )
}

# 분석할 지표 목록
metrics <- c("shannon", "chao1", "simpson", "faith_pd", "evenness", "observed_features")



# 각 지표에 대해 그래프 생성 및 저장
for (metric in metrics) {
    plot_and_save(metric, data, output_dir)
}
