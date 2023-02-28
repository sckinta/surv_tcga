# run_coxph_2d
run_coxph_2d <- function(formula, input_data, q1, q2, q1_col, q2_col){
    tmp <- input_data |> 
        select_at(all.vars(formula)) |> 
        na.omit() |> 
        mutate_at(q1_col, ~.x <= quantile(.x, q1)) |> 
        mutate_at(q2_col, ~.x <= quantile(.x, q2))
    
    tmp <- tmp |> 
        group_by_at(c(q1_col, q2_col)) |> 
        mutate(group=cur_group_id()) |> 
        ungroup() |> 
        mutate(group_label = case_when(
            group==1 ~ "high;high",
            group==2 ~ "low;high",
            group==3 ~ "high;low",
            group==4 ~ "low;low",
        )) |> 
        mutate(group=as.factor(group)) |> 
        I()
    
    
    
    x_sides <- attr(terms(formula), 'term.labels') # formula x side
    x_sides <- x_sides[!x_sides %in% c(q1_col, q2_col)]
    y_sides <- rownames(attr(terms(formula), 'factors'))[1] # forumula y side
    
    
    test <- coxph(as.formula(paste(y_sides, "~", "group +", x_sides)), data = tmp)
    
    list(test = test, data = tmp)
    
}

# optimal_cutoff_2d
optimal_cutoff_2d <- function(formula, input_data, q1_col, q2_col, low_q=0.2, high_q=0.8, breaks=10, pvalue_type=c("log", "sc", "wald")[1], min_sample_n=10){
    
    df <- pmap_dfr(
        crossing(
            q1 = seq(low_q, high_q, length.out=breaks),
            q2 = seq(low_q, high_q, length.out=breaks)
        ),
        function(q1, q2){
            res <- run_coxph_2d(formula, input_data, q1, q2, q1_col, q2_col)
            
            min_n <- res$data |> 
                group_by(group) |> 
                dplyr::count() |> 
                pull(n) |> min()
            
            res$test |> 
                broom::glance() |> 
                mutate(q1=q1, q2=q2) |> 
                mutate(cutoff1=quantile(expr, q1), cutoff2=quantile(expr, q2), min_sample_n=min_n)
        }
    ) |> 
        arrange(across(paste0("p.value.",pvalue_type))) |> 
        select_at(c("q1","cutoff1", "q2","cutoff2", paste0("p.value.",pvalue_type), "min_sample_n")) |> 
        mutate(keep = min_sample_n >= !!min_sample_n)
    
    q <- unlist(df[df$keep,][1,c("q1", "q2")])
    
    res <- run_coxph_2d(formula, input_data, q["q1"], q["q2"], q1_col, q2_col)
    
    list(q=unlist(df[df$keep,][1,c("q1", "q2")]), q_rank = df, best_coxph=res$test, best_data=res$data)
}


# plot_surv

plot_surv <-
    function(formula,
             data,
             strata_title = NULL,
             x_title = NULL,
             y_title = NULL,
             plot_title = NULL) {
        if (is.null(strata_title)) {
            strata_title = "strata"
        }
        if (is.null(x_title)) {
            x_title = "time"
        }
        if (is.null(y_title)) {
            y_title = "overall survival rate"
        }
        
        surv_fit2 <- survfit(formula,
                             data = data) # survfit obj
        
        surv_fit2 |>
            broom::tidy() |>
            left_join(
                broom::tidy(surv_fit2) |>
                    group_by(strata) |>
                    dplyr::slice(which.min(time)) |>
                    select(n = n.risk, strata)
            ) |>
            mutate(strata = str_replace(strata, ".*=", "")) |>
            filter(!grepl("\\[", strata)) |>
            mutate(strata = glue::glue("{strata} (n={n})")) |>
            ggplot(aes(x = time, y = estimate, group = strata)) +
            geom_line(aes(color = strata)) +
            geom_point(aes(color = strata), shape = "|", size = 2) +
            geom_ribbon(aes(
                ymax = conf.high,
                ymin = conf.low,
                fill = strata
            ), alpha = 0.2) +
            labs(x = x_title, y = y_title, color = strata_title) +
            guides(fill = "none") +
            ggtitle(plot_title)
    }

# optimal_cutoff_1d

optimal_cutoff_1d <- function(surv_obj, expr, low_q=0.2, high_q=0.8, breaks=10, log=T, pseudocount=1, pvalue_type=c("log", "sc", "wald")[1]){
    if(log){
        expr <- log2(pseudocount + expr)
    }
    
    df <- map_dfr(
        seq(low_q, high_q, length.out=breaks),
        ~coxph(surv_obj ~ (expr > quantile(expr, .x))) |> 
            broom::glance() |> 
            mutate(q=.x) |> 
            mutate(cutoff=quantile(expr, q))
    ) |> 
        arrange(across(paste0("p.value.",pvalue_type))) |> 
        select_at(c("q","cutoff", paste0("p.value.",pvalue_type)))
    df
}