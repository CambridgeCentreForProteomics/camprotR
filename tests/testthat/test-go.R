test_that("expand_terms() works", {
  # Imagine some made-up protein X which has the following GO terms annotated
  # mRNA cap binding, RNA binding, cytoplasmic region, innate immune response
  go_df <- data.frame(
    GO.ID =  c("GO:0098808", "GO:0003723", "GO:0099568", "GO:0045087"),
    GO.ONT = c("MF", "MF", "CC", "BP")
  )

  # need to run `get_all_mappings()` first
  expanded_go <- get_all_mappings(
    go_ids = go_df$GO.ID,
    ontologies = go_df$GO.ONT,
    verbose = TRUE,
    direction = "ancestor"
  )

  # construct input for `get_all_mappings()`
  go_2_ancestor <- go_df$GO.ID
  names(go_2_ancestor) <- go_df$GO.ONT

  # run `expand_terms()` to test and capture output
  # check output against snap reference
  expect_snapshot(
    expand_terms(
      go_df = go_df,
      go_col = "GO.ID",
      go2Ancestor = go_2_ancestor
    )
  )
})

## Generate some data to test other GO term functions
# get some significantly enriched proteins
dep_filt <- dep_results %>%
  distinct(ensembl, .keep_all = TRUE)

enriched <- dep_filt %>%
  mutate(sig_inc = Ubi4_vs_Ctrl_ratio > 0 & Ubi4_vs_Ctrl_p.adj < 0.05) %>%
  dplyr::select(ensembl, sig_inc) %>%
  tibble::deframe()

# make some fake bias data for generating pwf
set.seed(1234)
bias <- sample(1:25, size = length(enriched), replace = TRUE)

# generate pwf
pwf <- goseq::nullp(enriched, bias.data = bias, plot.fit = FALSE)

test_that("get_enriched_go() works", {
  # use dep_gene2cat which is installed with camprotR
  out <- get_enriched_go(
    pwf = pwf,
    gene2cat = dep_gene2cat
  ) %>%
    select(category, term)

  # check output against snap reference
  withr::with_tempfile("tf", {
    # save output to tempfile
    write.table(out, file = tf,
                sep = "\t", row.names = FALSE, col.names = TRUE)

    expect_snapshot_file(tf, "get-enriched-go.txt")
  })
})

test_that("estimate_go_overrep() works", {
  go_res <- get_enriched_go(
    pwf = pwf,
    gene2cat = dep_gene2cat
  )

  # does it work when gene2cat is named list?
  out <- estimate_go_overrep(
    obj = go_res,
    pwf = pwf,
    gene2cat = dep_gene2cat
  ) %>%
    select(category, term, adj_overrep) %>%
    mutate(adj_overrep = round(adj_overrep, digits = 3))

  # check output against snap reference
  withr::with_tempfile("tf", {
    # save output to tempfile
    write.table(out, file = tf,
                sep = "\t", row.names = FALSE, col.names = TRUE)

    expect_snapshot_file(tf, "estimate-go-overrep-list-input.txt")
  })

  gene_2_cat_df <- dep_gene2cat %>%
    tibble::enframe(name = "id", value = "go_term") %>%
    tidyr::unnest(cols = c("id", "go_term")) %>%
    as.data.frame()

  # does it work when gene2cat is a data.frame?
  out2 <- estimate_go_overrep(
    obj = go_res,
    pwf = pwf,
    gene2cat = gene_2_cat_df
  ) %>%
    select(category, term, adj_overrep) %>%
    mutate(adj_overrep = round(adj_overrep, digits = 3))

  # check output against snap reference
  withr::with_tempfile("tf", {
    # save output to tempfile
    write.table(out2, file = tf,
                sep = "\t", row.names = FALSE, col.names = TRUE)

    expect_snapshot_file(tf, "estimate-go-overrep-df-input.txt")
  })
})

test_that("remove_redundant_go() works", {
  go_res <- get_enriched_go(
    pwf = pwf,
    gene2cat = dep_gene2cat
  )

  # does it work when gene2cat is named list?
  go_res <- estimate_go_overrep(
    obj = go_res,
    pwf = pwf,
    gene2cat = dep_gene2cat
  )

  # check output against snap reference
  expect_snapshot(
    go_res %>%
      filter(grepl("immune", term)) %>%
      remove_redundant_go() %>%
      select(category, term, adj_overrep) %>%
      mutate(adj_overrep = round(adj_overrep, digits = 3))
  )
})

test_that("plot_go() works", {
  df <- data.frame(
    "term_short" = c("A GO term", "Another GO term", "One more GO term"),
    "ontology" = c("BP", "MF", "CC"),
    "over_represented_adj_pval" = c(0.0001, 1, 0.01),
    "adj_overrep" = c(15, 3, 1),
    "numDEInCat" = c(304, 22, 78)
  )

  vdiffr::expect_doppelganger(
    "plot-go",
    plot_go(df)
  )
})
