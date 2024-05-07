pacman::p_load(apila, tidymodels, here, tictoc, vroom, finetune, stacks, qs, glue, tidytable, tidyverse, lubridate, SMARTboost, xgboost, SHAPforxgboost, parallel, foreach, readr)

data_folder <- 'Set80'
save_folder <- 'Preds80_SMART'

data_path <- file.path(get_apila_repo_path(), 'Slices', data_folder)
save_path <- file.path(get_apila_repo_path(), 'Output', save_folder)

slices <- data_path %>%
  list.files() %>%
  keep(stringr::str_detect(., "Slice"))
last_slice <- slices[length(slices)] %>% stringr::str_remove(., 'model-') %>% stringr::str_remove(., '.rds')
#slices <- slices[c(7:11)]

splits <- qread(file.path(data_path, slices[1]), nthreads = 5)
test <- splits %>% testing %>% ungroup %>% select(-any_of('row_id'))

targets <-  test %>%
  select(contains('target')) %>%
  colnames() %>%
  enframe(value = 'target', name = 'purge') %>%
  mutate(horizon = sub("^[^_]*_", "", target)) %>%
  mutate(purge = case_when(stringr::str_detect(target, '1_2|3_4|2_3') ~ 1,
                           stringr::str_detect(target, '1_3|2_4|4_6') ~ 2,
                           stringr::str_detect(target, '1_4') ~ 3,
                           stringr::str_detect(target, '7_12') ~ 5,
                           TRUE ~ 0)) %>%
  expand_grid(fraction = c(0.67, 3/3),
              pval = c(0.1, 1)) %>%
  mutate(var = as.character(glue('pred_{horizon}kk_len{fraction}_feat{pval}'))) %>%
  filter(!(purge > 2 & fraction < 1))

rm(test, splits)

#haetaan kansiosta jo kertaalleen tuunatut parametrit
ready_tuned_params <- list.files(file.path(path_to_ml, 'Output', 'xgb_params'), full.names = FALSE) %>%
  enframe(value = 'var') %>%
  mutate(path = file.path(path_to_ml, 'Output', 'xgb_params', var),
         params = map(path, read_rds)) %>%
  select(var, params) %>%
  mutate(var = str_remove(var, '.rds'))

targets <- targets %>%
  #joinataan tuunatut parametrit
  left_join(ready_tuned_params) %>%
  fill(params, .direction = 'down')
targets

valid_metric <- 'ccc'
metrics <- metric_set(!!sym(valid_metric))
pimp_feat_selection <- FALSE #jos TRUE niin ottaa kaiken
train_len <- TRUE #jos TRUE niin ottaa kaiken
era_boost <- FALSE
neutralize_features <- FALSE
top <- c(0.1) #neutralization tops
prop <- c(0.5) #neutralization props


if(!pimp_feat_selection){
  targets <- targets %>% filter(pval == 1)
}
if(!train_len){
  targets <- targets %>% filter(fraction < 1)
}
treenattavat <- c('target_1_scv')

targets <- targets %>%
  filter(target %in% treenattavat)
#valiaikainen
# targets <- targets %>% dplyr::slice(c(1,4,9,14))

# targets <- targets %>%
#   bind_rows(targets %>% dplyr::slice(1) %>% mutate(target = 'target_1_scvol')) %>%
#   mutate(horizon = sub("^[^_]*_", "", target)) %>%
#   mutate(var = as.character(glue('pred_{horizon}kk_len{fraction}_feat{pval}')))

options(tidymodels.dark = TRUE)

for (slice in slices) {

  tryCatch(
    {
      gc()

      message("Reading file...")
      id <- slice %>% stringr::str_remove(., 'model-') %>% stringr::str_remove(., '.rds')

      splits <- qread(file.path(data_path, glue::glue('{slice}')), nthreads = 5)
      main_train <- splits %>% training %>% ungroup %>% select(-any_of('row_id'))
      test <- splits %>% testing %>% ungroup %>% select(-any_of('row_id'))

      pryr::mem_change(rm(splits))

      for (i in seq_len(nrow(targets))) {

        target_var <- targets$target[i]
        n_purge <- targets$purge[i]
        horizon <- targets$horizon[i]
        train_frac <- targets$fraction[i]
        pval <- targets$pval[i]
        var <- sym(glue('pred_{horizon}kk_len{train_frac}_feat{pval}'))

        date_filter <- main_train %>%
          distinct(date) %>%
          mutate(row = row_number()) %>%
          filter(row >= quantile(row, 1-train_frac)) %>%
          pull(date) %>%
          first()

        train <- exclude_targets_from_data(target_col = target_var, data = filter(main_train, date >= date_filter)) %>%
          filter_tails(target = target_var, tail = 0.2) %>%
          purge_dataset(n = n_purge) %>%
          ungroup

        #Joka targetille feat selection
        if(pval < 1){

          message(glue('Selecting features'))
          features <- shap_imp(train, target_var, S = 50, parallel = TRUE, sample = TRUE, sample_n = 0.2,
                               ncores = 2, seed = 123, tree_method = "gpu_hist")

          selected_features <- features %>%
            filter(p_value <= pval) %>%
            pull(var)
          n_features <- length(selected_features)

          if(n_features > 0){
            train <- train %>%
              select(Id, date, contains('target'), any_of(selected_features))
          }

          rm(features, selected_features)
        } else {
          n_features <- ncol(train)
        }

        # message(glue('{id}, training {target_var}, train length {train_frac} with {n_features} features'))

        # resamples <- train %>%
        #   arrange(date) %>%
        #   generate_ts_resamples(max_samples = 6, period = 'month') %>%
        #   validation_purge(n_purge)
        # gc()
        #
        train <- train %>%
          mutate(across(where(is.numeric) & !starts_with("target"), ~./10000))

        # dates <- train %>%
        #   select(date) %>%
        #   distinct() %>%
        #   pull(date)
        #
        # y <- train[,"target_1"] - mean(train[,"target_1"]$target_1)
        #
        # ssc <- 0.0
        #
        # for(date in dates){
        #   ssc = ssc + (sum(y[train$date == date,]))^2
        # }
        #
        # loglikdivide  = ssc/sum(y$target_1^2)
        #
        SMART_recipe <- recipe(train) %>%
          update_role(everything() & !any_of('case_wts')) %>%
          update_role(all_of(target_var), new_role = "outcome") %>%
          #step_impute_median(all_numeric_predictors()) %>%
          update_role(any_of('Id'), new_role = "ID") %>%
          step_unknown(all_nominal_predictors(), new_level='missing') %>% # replace missing characters with unknown
          step_novel(all_nominal_predictors(), new_level='Unseen') %>% # deal with new data being present in new data
          step_other(all_nominal_predictors(), other='Misc') %>% # collapse infrequent levels
          step_rm(contains("row") | contains("macro") | contains("fct") | contains("sampl_")) %>% #otetaan n?m? pois
          step_nzv(all_numeric_predictors()) %>%
          step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
          step_impute_median(all_numeric_predictors())

        SMART_spec  <- SMARTboost(
          trees = 50,
          mtry = 0.6,
          tree_depth = 4,
          learn_rate = 0.2)

        SMART_spec <- SMART_spec %>%
          set_engine("SMARTboost", optimizevs = FALSE) %>%
          set_mode("regression")


        wf <- workflow() %>%
          add_recipe(SMART_recipe) %>%
          add_model(SMART_spec)

        fit <- fit(wf, train)

        suppressMessages(
          target_preds <- test %>%
            group_by(date) %>%
            group_split() %>%
            map_dfr(~predict(fit, .)) %>%
            transmute("{var}" := rowMeans(across(contains('pred')), na.rm = TRUE))
        )
        test <- test %>%
          bind_cols(target_preds)

        pryr::mem_change(rm(wf, resamples, train, target_preds))


        test <- test %>%
          select(date, Id, contains('pred'), everything()) %>%
          qsave(file = file.path(save_path, paste0(id, '.rds')))

        pryr::mem_change(rm(main_train))
      }

    },
error=function(error_message) {
  message(error_message)
}
)
}
