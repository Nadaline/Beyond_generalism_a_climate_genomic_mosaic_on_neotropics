###############################################################################
# GRADIENT FOREST COM TODOS OS 14.134 SNPs
# Ameiva ameiva: modelo bruto e modelo controlado por estrutura + geografia
#
# COMO RODAR NO TERMINAL (a partir da pasta do projeto):
# nohup Rscript GF_14134_bruto_controlado_8cores.R > GF_14134_live.log 2>&1 &
# tail -f GF_14134_live.log
#
# O script usa ate 8 nucleos, salva cada bloco no HD externo e pode continuar
# de onde parou. Ele NAO altera os resultados antigos do RDA-2027.R.
###############################################################################


# 1. PACOTES ####

# Inclui primeiro a biblioteca local do projeto, onde gradientForest esta.
.libPaths(c(normalizePath("R-library"), .libPaths()))

# gradientForest ajusta uma Random Forest separada para cada SNP.
library(gradientForest)

# parallel permite que diferentes blocos de SNPs sejam rodados simultaneamente.
library(parallel)

# LEA le a matriz .lfmm imputada que ja foi produzida no RDA-2027.
library(LEA)


# 2. DECISOES QUE VOCE PODE ALTERAR ####

# Numero maximo de processos simultaneos. Em seu computador, usar 8.
N_CORES <- 6L

# Numero de arvores por SNP. 500 e apropriado para a primeira execucao completa.
# Depois, os SNPs mais estaveis podem ser reajustados com 1.000 ou 2.000.
N_TREES <- 500L

# Numero de variaveis ambientais candidatas em cada divisao da arvore.
# Com dez preditores, mtry = 2 oferece um meio-termo adequado.
MTRY <- 3L

# Usamos mais blocos do que nucleos. Assim, cada arquivo e menor e, se a
# execucao cair, perde-se no maximo um bloco curto, nao 1/8 de toda a analise.
N_BLOCKS <- 32L

# Semente geral para que a divisao dos SNPs e os forests sejam reproduziveis.
MASTER_SEED <- 2027L

# TRUE roda o GF com genotipos 0/1/2 sem controle de estrutura.
RUN_RAW <- TRUE

# TRUE roda o GF sobre residuos controlados por PC1-PC3, longitude e latitude.
RUN_CONTROLLED <- TRUE

# TRUE combina oficialmente os blocos com combinedGradientForest().
COMBINE_BLOCKS <- TRUE

# Se TRUE, arquivos de blocos ja completos sao reutilizados.
RESUME_FINISHED_BLOCKS <- TRUE


# 3. CAMINHOS ####

project <- "/home/nadaline/Documents/Phd/RADseq/The-Edge-of-Adaptation"

# Os dados grandes e todos os novos resultados ficam no HD externo.
base_results <- paste0(
  "/media/nadaline/NADALINE II/The-Edge-of-Adaptation/",
  "genome_annotation/results/RDA-2027-14134"
)

# Pasta exclusiva desta analise, para nao sobrescrever o GF de 118 candidatos.
GF_out <- file.path(base_results, "gradient_forest_all_14134")
GF_raw_out <- file.path(GF_out, "raw")
GF_controlled_out <- file.path(GF_out, "controlled_PC1_PC3_lon_lat")

dir.create(GF_raw_out, recursive = TRUE, showWarnings = FALSE)
dir.create(GF_controlled_out, recursive = TRUE, showWarnings = FALSE)

# Arquivos de entrada ja gerados e auditados pelo RDA-2027.
genotype_cache <- file.path(base_results, "03_genotype_matrix_14134.rds")
imputed_file <- file.path(base_results, "03_Ameiva_14134.lfmm_imputed.lfmm")
pca_file <- file.path(base_results, "04_PCA_14134_scores.csv")
environment_decision_file <- file.path(
  base_results, "06_environment_forward_selection.csv"
)
environment_file <- "/home/nadaline/Documents/Phd/RADseq/Ameiva_Variables_rda.csv"

# Interrompe cedo, com mensagem clara, caso o HD nao esteja montado.
required_files <- c(
  genotype_cache, imputed_file, pca_file,
  environment_decision_file, environment_file
)
stopifnot(all(file.exists(required_files)))

# Nao permite solicitar mais nucleos do que o computador oferece.
available_cores <- parallel::detectCores(logical = FALSE)
N_CORES <- min(N_CORES, available_cores, N_BLOCKS)


# 4. MATRIZ GENOTIPICA IMPUTADA: 245 INDIVIDUOS x 14.134 SNPs ####

# O cache contem nomes dos individuos, IDs dos SNPs e coordenadas genomicas.
genotype_object <- readRDS(genotype_cache)

# A matriz .lfmm contem os mesmos genotipos depois da imputacao por sNMF K=4.
geno_complete <- as.matrix(LEA::read.lfmm(imputed_file))

# Recoloca nomes de linhas e colunas, pois o formato LFMM nao os armazena.
rownames(geno_complete) <- rownames(genotype_object$geno)
colnames(geno_complete) <- colnames(genotype_object$geno)

# Auditorias essenciais: o GF nao deve iniciar com painel ou ordem incorretos.
stopifnot(nrow(geno_complete) == 245L)
stopifnot(ncol(geno_complete) == 14134L)
stopifnot(!anyNA(geno_complete))
stopifnot(all(geno_complete %in% c(0, 1, 2)))


# 5. AMBIENTE E ALINHAMENTO DOS INDIVIDUOS ####

# check.names=FALSE preserva os nomes originais CHELSA.
Variables <- read.csv(environment_file, check.names = FALSE)

# Usa InternalID como chave e reorganiza o CSV na ordem exata dos genotipos.
rownames(Variables) <- Variables$InternalID
Variables <- Variables[rownames(geno_complete), , drop = FALSE]
stopifnot(identical(rownames(Variables), rownames(geno_complete)))

# Recupera do próprio RDA-2027 quais preditores foram selecionados.
#environment_decision <- read.csv(environment_decision_file)
#selected_environment <- strsplit(
#  environment_decision$selected[1], ";", fixed = TRUE
#)[[1]]
#selected_environment <- trimws(selected_environment)
# Uma recomendação para America do sul e não usarmos as variaveis trimestrais que tem problemas de continuidade espacial.
# Lista de bioclimas recomendada (sem trimestres e sem Bio14)
selected_environment <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio12", "bio13", "bio15")

# Converte nomes curtos como bio2 para as colunas CHELSA do CSV.
environment_columns <- paste0( "X.CHELSA_", selected_environment, "_1981.2010_SA")
stopifnot(all(environment_columns %in% names(Variables)))

# Copia somente as variaveis escolhidas e restaura os nomes biologicos curtos.
Environment_raw <- Variables[, environment_columns, drop = FALSE]
names(Environment_raw) <- selected_environment
stopifnot(!anyNA(Environment_raw))

# ECDF coloca cada preditor entre 0 e 1, preserva a ordem dos valores e evita
# que poucos extremos controlem as divisoes das arvores.
Environment <- as.data.frame(lapply(
  Environment_raw,
  function(x) stats::ecdf(x)(x)
))
rownames(Environment) <- rownames(geno_complete)


# 6. COVARIAVEIS DE ESTRUTURA E GEOGRAFIA ####

# Carrega os PCs calculados exatamente sobre o painel de 14.134 SNPs.
PCA_scores <- read.csv(pca_file, check.names = FALSE)
rownames(PCA_scores) <- PCA_scores$InternalID
PCA_scores <- PCA_scores[rownames(geno_complete), , drop = FALSE]

# Longitude e latitude sao padronizadas para terem escalas comparaveis.
covariates <- data.frame(
  PC1 = PCA_scores$PC1,
  PC2 = PCA_scores$PC2,
  PC3 = PCA_scores$PC3,
  longitude = as.numeric(scale(Variables$long)),
  latitude = as.numeric(scale(Variables$lat))
)
rownames(covariates) <- rownames(geno_complete)
stopifnot(!anyNA(covariates))

# A matriz inclui intercepto. QR e calculada uma unica vez para todos os SNPs.
covariate_matrix <- model.matrix(
  ~ PC1 + PC2 + PC3 + longitude + latitude,
  data = covariates
)
covariate_qr <- qr(covariate_matrix)

# Remove de cada SNP a parte linear explicada por estrutura e geografia.
# qr.resid aceita toda a matriz de SNPs de uma vez, sendo muito mais rapido do
# que chamar lm() 14.134 vezes.
geno_controlled <- qr.resid(covariate_qr, geno_complete)
rownames(geno_controlled) <- rownames(geno_complete)
colnames(geno_controlled) <- colnames(geno_complete)

# Os residuos sao centrados; padronizamos apenas a variancia para que SNPs de
# frequencias diferentes nao tenham peso determinado apenas por sua escala.
residual_sd <- apply(geno_controlled, 2, sd)
stopifnot(all(is.finite(residual_sd)), all(residual_sd > 0))
geno_controlled <- sweep(geno_controlled, 2, residual_sd, "/")


# 7. MAPA ENTRE NOMES INTERNOS DO GF E IDs REAIS DOS SNPs ####

# Nomes simples evitam erros do pacote com IDs contendo dois-pontos e sinais.
GF_names <- sprintf("SNP%05d", seq_len(ncol(geno_complete)))
GF_map <- data.frame(
  GF_name = GF_names,
  stacks_id = colnames(geno_complete),
  scaffold = genotype_object$variant_table$CHROM,
  position = genotype_object$variant_table$POS,
  rad_locus = sub(":.*$", "", colnames(geno_complete)),
  stringsAsFactors = FALSE
)
write.csv(GF_map, file.path(GF_out, "00_GF_SNP_map.csv"), row.names = FALSE)

# As duas matrizes recebem os mesmos nomes e, portanto, podem ser comparadas.
colnames(geno_complete) <- GF_names
colnames(geno_controlled) <- GF_names


# 8. DIVISAO REPRODUZIVEL EM BLOCOS ####

# round-robin deixa todos os blocos praticamente do mesmo tamanho.
block_id <- rep(seq_len(N_BLOCKS), length.out = length(GF_names))
GF_blocks <- split(GF_names, block_id)

write.csv(
  data.frame(
    GF_name = GF_names,
    block = block_id,
    stacks_id = GF_map$stacks_id
  ),
  file.path(GF_out, "00_GF_block_registry.csv"),
  row.names = FALSE
)


# 9. FUNCAO QUE AJUSTA UM BLOCO ####

fit_GF_block <- function(block_number, response_matrix, model_name, model_out) {

  # Nome deterministico permite retomar a execucao automaticamente.
  block_file <- file.path(
    model_out, sprintf("GF_%s_block_%02d_of_%02d.rds",
                      model_name, block_number, N_BLOCKS)
  )

  # Se o bloco ja existe e pode ser lido, nao gastamos horas recalculando-o.
  if (RESUME_FINISHED_BLOCKS && file.exists(block_file)) {
    previous <- try(readRDS(block_file), silent = TRUE)
    if (!inherits(previous, "try-error") &&
        inherits(previous, "gradientForest")) {
      message("Reutilizando bloco pronto: ", basename(block_file))
      return(block_file)
    }
  }

  block_SNPs <- GF_blocks[[block_number]]

  # Cada bloco contem todos os individuos, os quatro preditores e parte dos SNPs.
  block_data <- data.frame(
    Environment,
    response_matrix[, block_SNPs, drop = FALSE],
    check.names = FALSE
  )

  # Cada bloco recebe semente propria, mas reproduzivel.
  set.seed(MASTER_SEED + block_number +
             ifelse(model_name == "controlled", 10000L, 0L))

  message(
    "Iniciando ", model_name, " | bloco ", block_number, "/", N_BLOCKS,
    " | ", length(block_SNPs), " SNPs | PID ", Sys.getpid()
  )

  block_model <- gradientForest(
    data = block_data,
    predictor.vars = names(Environment),
    response.vars = block_SNPs,
    ntree = N_TREES,
    mtry = MTRY,
    maxLevel = 0,
    # Com maxLevel=0, corr.threshold nao filtra respostas nem preditores.
    corr.threshold = .50,
    compact = TRUE,
    nbin = 101,
    trace = FALSE,
    check.names = FALSE
  )

  # Primeiro grava em arquivo temporario e depois renomeia. Isso impede que uma
  # interrupcao deixe um .rds incompleto parecendo um bloco valido.
  temporary_file <- paste0(block_file, ".partial")
  saveRDS(block_model, temporary_file, compress = FALSE)
  file.rename(temporary_file, block_file)

  message("Finalizado: ", basename(block_file))
  block_file
}


# 10. FUNCAO QUE RODA OS 32 BLOCOS USANDO ATE 8 NUCLEOS ####

run_GF_model <- function(response_matrix, model_name, model_out) {

  message("\n============================================================")
  message("MODELO: ", toupper(model_name))
  message("Nucleos: ", N_CORES, " | blocos: ", N_BLOCKS,
          " | arvores/SNP: ", N_TREES)
  message("Inicio: ", Sys.time())
  message("============================================================\n")

  # mclapply usa processos fork no Linux. No maximo oito blocos rodam juntos.
  block_files <- parallel::mclapply(
    seq_len(N_BLOCKS),
    fit_GF_block,
    response_matrix = response_matrix,
    model_name = model_name,
    model_out = model_out,
    mc.cores = N_CORES,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
  block_files <- unlist(block_files, use.names = FALSE)

  # Confirma que todos os blocos existem antes de resumir ou combinar.
  stopifnot(length(block_files) == N_BLOCKS)
  stopifnot(all(file.exists(block_files)))

  block_models <- lapply(block_files, readRDS)
  stopifnot(all(vapply(block_models, inherits, logical(1), "gradientForest")))

  # Extrai R2 fora-da-bolsa de cada SNP. O pacote mantem em result somente os
  # SNPs com R2 positivo; os demais recebem zero na tabela completa.
  positive_R2 <- unlist(lapply(block_models, function(x) x$result))
  positive_names <- sub("^[^.]+\\.", "", names(positive_R2))
  if (length(positive_names) == length(positive_R2)) names(positive_R2) <- positive_names

  SNP_summary <- GF_map
  SNP_summary$R2_OOB <- 0
  matched <- match(names(positive_R2), SNP_summary$GF_name)
  valid <- !is.na(matched)
  SNP_summary$R2_OOB[matched[valid]] <- as.numeric(positive_R2[valid])
  SNP_summary$positive_R2 <- SNP_summary$R2_OOB > 0
  SNP_summary <- SNP_summary[order(SNP_summary$R2_OOB, decreasing = TRUE), ]

  write.csv(
    SNP_summary,
    file.path(model_out, paste0("GF_", model_name, "_SNP_R2.csv")),
    row.names = FALSE
  )

  # Importancia ambiental media ponderada igualmente entre os blocos.
  block_importance <- do.call(rbind, lapply(
    seq_along(block_models),
    function(i) data.frame(
      block = i,
      predictor = names(block_models[[i]]$overall.imp),
      importance = as.numeric(block_models[[i]]$overall.imp)
    )
  ))
  write.csv(
    block_importance,
    file.path(model_out, paste0("GF_", model_name,
                               "_predictor_importance_by_block.csv")),
    row.names = FALSE
  )

  # A funcao oficial combina curvas cumulativas e permite predict() posterior.
  combined_model <- NULL
  if (COMBINE_BLOCKS) {
    names(block_models) <- sprintf("block%02d", seq_along(block_models))
    combined_model <- do.call(
      gradientForest::combinedGradientForest,
      c(block_models, list(nbin = 101, method = 2, standardize = "before"))
    )
    saveRDS(
      combined_model,
      file.path(model_out, paste0("GF_", model_name,
                                 "_combined_14134.rds")),
      compress = FALSE
    )
  }

  # Libera os modelos individuais da memoria; os arquivos permanecem no HD.
  rm(block_models)
  gc()

  message("Termino de ", model_name, ": ", Sys.time())
  list(SNP_summary = SNP_summary, combined_model = combined_model)
}


# 11. EXECUCAO BRUTA E CONTROLADA ####

set.seed(MASTER_SEED)

GF_raw <- NULL
if (RUN_RAW) {
  GF_raw <- run_GF_model(
    response_matrix = geno_complete,
    model_name = "raw",
    model_out = GF_raw_out
  )
}

GF_controlled <- NULL
if (RUN_CONTROLLED) {
  GF_controlled <- run_GF_model(
    response_matrix = geno_controlled,
    model_name = "controlled",
    model_out = GF_controlled_out
  )
}


# 12. COMPARACAO ENTRE GF BRUTO, GF CONTROLADO, RDA E LFMM ####


if (RUN_RAW && RUN_CONTROLLED) {
  
  comparison <- merge(
    GF_raw$SNP_summary[, c(
      "GF_name", "stacks_id", "rad_locus", "R2_OOB"
    )],
    GF_controlled$SNP_summary[, c("GF_name", "R2_OOB")],
    by = "GF_name",
    suffixes = c("_raw", "_controlled"),
    all = TRUE
  )
  
  # A tabela deve conter o painel inteiro, não apenas candidatos prévios.
  stopifnot(nrow(comparison) == nrow(GF_map))
  stopifnot(!anyDuplicated(comparison$GF_name))
  
RDA_file <- file.path(
  base_results,
  "08_RDA_candidates_landscape_consensus.csv"
)

RDA_candidates <- read.csv(RDA_file, check.names = FALSE)$stacks_id
comparison$RDA <- comparison$stacks_id %in% RDA_candidates
  
  RDA_candidates <- read.csv(RDA_file, check.names = FALSE)$stacks_id
  comparison$RDA <- comparison$stacks_id %in% RDA_candidates
  # Apenas o resultado LFMM2 atual entra na integracao. O arquivo LFMM antigo
  # (ridge sobre PCs ambientais) e mantido somente como registro historico.
  LFMM_file <- file.path(
    base_results, "09_LFMM2_candidates_FDR05.csv"
  )
  

  
  LFMM_candidates <- if (file.exists(LFMM_file)) {
    read.csv(LFMM_file, check.names = FALSE)$stacks_id
  } else {
    character()
  }
  
  comparison$RDA <- comparison$stacks_id %in% RDA_candidates
  comparison$LFMM <- comparison$stacks_id %in% LFMM_candidates
  
  # “R2 positivo” não é um filtro útil aqui. Definimos alto suporte GF
  # de forma transparente: 5% superiores do GF controlado.
  GF_R2_threshold <- quantile(
    comparison$R2_OOB_controlled,
    probs = 0.95,
    na.rm = TRUE,
    names = FALSE
  )
  
  comparison$GF_high <- comparison$R2_OOB_controlled >= GF_R2_threshold
  comparison$linear_support <- comparison$RDA | comparison$LFMM
  
  # As classes descrevem padrões de detecção; não são um ranking de qualidade.
  comparison$detection_class <- ifelse(
    comparison$GF_high & !comparison$linear_support,
    "GF prioritário",
    ifelse(
      !comparison$GF_high & comparison$linear_support,
      "Linear prioritário",
      ifelse(
        comparison$GF_high & comparison$linear_support,
        "Suporte compartilhado",
        "Sem suporte prioritário"
      )
    )
  )
  
  comparison$detection_class <- factor(
    comparison$detection_class,
    levels = c(
      "GF prioritário",
      "Linear prioritário",
      "Suporte compartilhado",
      "Sem suporte prioritário"
    )
  )
  
  write.csv(
    comparison,
    file.path(GF_out, "GF_RDA_LFMM_full_14134_comparison.csv"),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(
      metric = c("total_SNPs", "GF_R2_percentile", "GF_R2_threshold"),
      value = c(nrow(comparison), 0.95, GF_R2_threshold)
    ),
    file.path(GF_out, "GF_RDA_LFMM_classification_settings.csv"),
    row.names = FALSE
  )
}


# 13. REGISTRO FINAL ####

run_summary <- data.frame(
  item = c(
    "date", "individuals", "SNPs", "RAD_loci", "environment",
    "covariates_controlled_GF", "cores", "blocks", "trees_per_SNP",
    "raw_completed", "controlled_completed"
  ),
  value = c(
    as.character(Sys.time()), nrow(geno_complete), ncol(geno_complete),
    length(unique(GF_map$rad_locus)), paste(names(Environment), collapse = ";"),
    "PC1;PC2;PC3;longitude;latitude", N_CORES, N_BLOCKS, N_TREES,
    RUN_RAW, RUN_CONTROLLED
  )
)
write.csv(run_summary, file.path(GF_out, "00_GF_run_summary.csv"), row.names = FALSE)

message("\nGradient Forest bruto e controlado finalizados com sucesso.")
message("Resultados: ", GF_out)
###############################################################################
##### Avaliation
base_results <- "/media/nadaline/NADALINE II/The-Edge-of-Adaptation/genome_annotation/results/RDA-2027-14134"
GF_out <- file.path(base_results, "gradient_forest_all_14134")

gf_raw <- readRDS(file.path(GF_out, "raw/GF_raw_combined_14134.rds"))
gf_controlled <- readRDS(file.path(
  GF_out,
  "controlled_PC1_PC3_lon_lat/GF_controlled_combined_14134.rds"
))
plot(gf_raw, plot.type = "Overall.Importance")
plot(gf_raw, plot.type = "Cumulative.Importance")

plot(gf_controlled, plot.type = "Overall.Importance")
plot(gf_controlled, plot.type = "Cumulative.Importance")

raw_r2 <- read.csv(file.path(GF_out, "raw/GF_raw_SNP_R2.csv"))

hist(raw_r2$R2_OOB,
     breaks = 60,
     main = "Distribuição de R² OOB — GF bruto",
     xlab = "R² OOB por SNP")

head(raw_r2[order(-raw_r2$R2_OOB), ], 30)

comparison <- read.csv(
  file.path(GF_out, "GF_raw_controlled_RDA_LFMM_SNP_comparison.csv")
)

strong <- subset(
  comparison,
  GF_controlled_positive & (RDA | LFMM)
)

strong[order(-strong$R2_OOB_controlled), ]


purple <- "#6B2D8B"
orange <- "#E87516"

pdf(
  file.path(GF_out, "Figure_GF_cumulative_raw_controlled.pdf"),
  width = 12,
  height = 6
)

par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1))

plot(
  gf_raw,
  plot.type = "Cumulative.Importance",
  col = purple,
  main = "A. GF bruto"
)

plot(
  gf_controlled,
  plot.type = "Cumulative.Importance",
  col = orange,
  main = "B. GF controlado"
)

(mtext(
  "Mudança genômica cumulativa ao longo dos gradientes climáticos",
  side = 3, outer = TRUE, line = -1.2
))

dev.off()



##### Final aproach

library(ggplot2)

COL_PURPLE <- "#6B3F92"
COL_ORANGE <- "#F39940"
COL_TEAL   <- "#297C6F"

plot_data <- subset(
  comparison,
  detection_class != "Sem suporte prioritário"
)

p_GF <- ggplot(
  plot_data,
  aes(
    x = R2_OOB_raw,
    y = R2_OOB_controlled,
    colour = detection_class
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linewidth = 0.7,
    linetype = 2,
    colour = "grey45"
  ) +
  geom_point(size = 1.6, alpha = 0.65) +
  coord_equal(xlim = c(0, 0.8), ylim = c(0, 0.8)) +
  scale_colour_manual(
    values = c(
      "GF prioritário" = COL_TEAL,
      "Linear prioritário" = COL_PURPLE,
      "Suporte compartilhado" = COL_ORANGE
    ),
    name = "Classe de detecção",
    labels = c(
      "GF prioritário" =
        "GF prioritário: hipótese de resposta não linear/condicional",
      "Linear prioritário" =
        "Linear prioritário: RDA e/ou LFMM",
      "Suporte compartilhado" =
        "Suporte compartilhado: GF + RDA e/ou LFMM"
    )
  ) +
  labs(
    x = expression(R^2~"OOB — GF bruto"),
    y = expression(R^2~"OOB — GF controlado")
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = c(0.68, 0.84),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

ggsave(
  file.path(GF_out, "Figure_GF_raw_controlled_detection_classes.png"),
  p_GF,
  width = 9,
  height = 7,
  dpi = 600
)


##### trÊs conjuntos de snps candidatos, lineares, nao lineares e sem o gf
#### CLASSIFICACAO DOS TRES CONJUNTOS PRIORITARIOS ####

# Limiar explícito para “alto suporte” no GF controlado.
GF_R2_threshold <- quantile(
  comparison$R2_OOB_controlled,
  probs = 0.95,
  na.rm = TRUE,
  names = FALSE
)

comparison$GF_high <- comparison$R2_OOB_controlled >= GF_R2_threshold

# Mantém todas as combinações na tabela, inclusive as intermediárias.
comparison$candidate_set <- "Intermediário / outro"

# 1. Consenso estrito: menor risco de falso positivo.
comparison$candidate_set[
  comparison$GF_high & comparison$RDA & comparison$LFMM
] <- "Consenso triplo"

# 2. GF prioritário: hipótese de relação não linear ou condicional.
comparison$candidate_set[
  comparison$GF_high & !comparison$RDA & !comparison$LFMM
] <- "GF prioritário"

# 3. Linear prioritário: efeito recuperado por RDA e/ou LFMM,
# mas sem alto suporte relativo no GF.
comparison$candidate_set[
  !comparison$GF_high & (comparison$RDA | comparison$LFMM)
] <- "Linear prioritário"

comparison$candidate_set <- factor(
  comparison$candidate_set,
  levels = c(
    "GF prioritário",
    "Linear prioritário",
    "Consenso triplo",
    "Intermediário / outro"
  )
)

# Conferencia dos tamanhos dos conjuntos.
candidate_set_summary <- as.data.frame(
  table(comparison$candidate_set),
  stringsAsFactors = FALSE
)
names(candidate_set_summary) <- c("candidate_set", "n_SNPs")

write.csv(
  comparison,
  file.path(GF_out, "GF_RDA_LFMM_full_14134_comparison.csv"),
  row.names = FALSE
)

write.csv(
  candidate_set_summary,
  file.path(GF_out, "GF_RDA_LFMM_candidate_set_summary.csv"),
  row.names = FALSE
)

# Arquivo prático apenas com os três conjuntos que serão interpretados.
priority_sets <- subset(
  comparison,
  candidate_set != "Intermediário / outro"
)

write.csv(
  priority_sets,
  file.path(GF_out, "GF_RDA_LFMM_three_priority_sets.csv"),
  row.names = FALSE
)

#ploting
COL_PURPLE <- "#6B3F92"
COL_ORANGE <- "#F39940"
COL_TEAL   <- "#297C6F"

plot_data <- subset(
  comparison,
  candidate_set %in% c(
    "GF prioritário",
    "Linear prioritário",
    "Consenso triplo"
  )
)

p_GF <- ggplot(
  plot_data,
  aes(
    x = R2_OOB_raw,
    y = R2_OOB_controlled,
    colour = candidate_set
  )
) +
  geom_abline(
    slope = 1, intercept = 0,
    linewidth = 0.7,
    linetype = 2,
    colour = "grey45"
  ) +
  geom_point(size = 1.6, alpha = 0.65) +
  coord_equal(xlim = c(0, 0.9), ylim = c(0, 0.9)) +
  scale_colour_manual(
    values = c(
      "GF prioritário" = COL_TEAL,
      "Linear prioritário" = COL_PURPLE,
      "Consenso triplo" = COL_ORANGE
    ),
    name = "Conjunto prioritário",
    labels = c(
      "GF prioritário" =
        "GF prioritário: hipótese não linear/condicional",
      "Linear prioritário" =
        "Linear prioritário: RDA e/ou LFMM",
      "Consenso triplo" =
        "Consenso triplo: GF + RDA + LFMM"
    )
  ) +
  labs(
    x = expression(R^2~"OOB — GF bruto"),
    y = expression(R^2~"OOB — GF controlado")
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

ggsave(
  file.path(GF_out, "Figure_GF_three_priority_sets.png"),
  p_GF,
  width = 9,
  height = 7,
  dpi = 600
)

#Acima basicamente vemos em laranja um conjunto de preditores lineares corroborado pelas 3 analises, como GF analise nao linearidade resulta em os loci nao compartilhados em verde, que pelo r² sao os candidatos mais promissores. Vamos ver isso em um vennplot:
#### VENN: GF CONTROLADO, RDA E LFMM ####

library(ggVennDiagram)

# O GF entra como o conjunto de loci do percentil 95 de R² controlado.
GF_R2_threshold <- quantile(
  comparison$R2_OOB_controlled,
  probs = 0.95,
  na.rm = TRUE,
  names = FALSE
)

GF_high_ids <- comparison$stacks_id[
  comparison$R2_OOB_controlled >= GF_R2_threshold
]

RDA_ids <- comparison$stacks_id[comparison$RDA]
LFMM_ids <- comparison$stacks_id[comparison$LFMM]

venn_sets <- list(
  "GF controlado\n(top 5% R² OOB)" = GF_high_ids,
  "RDA" = RDA_ids,
  "LFMM" = LFMM_ids
)

p_venn <- ggVennDiagram(
  venn_sets,
  label_alpha = 0
) +
  scale_fill_gradient(
    low = "#F8F5F2",
    high = "#6B3F92"
  ) +
  scale_color_manual(
    values = c(
      "GF controlado\n(top 5% R² OOB)" = "#297C6F",
      "RDA" = "#6B3F92",
      "LFMM" = "#F39940"
    )
  ) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Complementaridade entre abordagens de associação genótipo–ambiente"
  )

ggsave(
  file.path(GF_out, "Figure_Venn_GFcontrolled_RDA_LFMM.png"),
  p_venn,
  width = 8,
  height = 7,
  dpi = 600
)

p_venn

library(VennDiagram)
library(grid)

venn_plot <- draw.triple.venn(
  area1 = length(GF_high_ids),
  area2 = length(RDA_ids),
  area3 = length(LFMM_ids),
  n12 = length(intersect(GF_high_ids, RDA_ids)),
  n13 = length(intersect(GF_high_ids, LFMM_ids)),
  n23 = length(intersect(RDA_ids, LFMM_ids)),
  n123 = length(Reduce(intersect, venn_sets)),
  category = c("GF controlado\n(top 5% R² OOB)", "RDA", "LFMM"),
  fill = c("#297C6F", "#6B3F92", "#F39940"),
  alpha = c(0.40, 0.40, 0.40),
  lty = "blank",
  cex = 1.5,
  cat.cex = 1.3,
  cat.col = c("#297C6F", "#6B3F92", "#F39940")
)

png(
  file.path(GF_out, "Figure_Venn_GFcontrolled_RDA_LFMM_colored.png"),
  width = 2400,
  height = 2000,
  res = 300
)

grid.newpage()
grid.draw(venn_plot)
dev.off()
