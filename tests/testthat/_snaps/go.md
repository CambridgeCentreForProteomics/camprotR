# expand_terms() works

    Code
      expand_terms(go_df = go_df, go_col = "GO.ID", go2Ancestor = go_2_ancestor)
    Output
              GOID
      1 GO:0098808
      2 GO:0003723
      3 GO:0099568
      4 GO:0045087
      5       <NA>

# get_enriched_go() works

    Code
      lapply(out, class)
    Output
      $category
      [1] "character"
      
      $over_represented_pvalue
      [1] "numeric"
      
      $under_represented_pvalue
      [1] "numeric"
      
      $numDEInCat
      [1] "integer"
      
      $numInCat
      [1] "integer"
      
      $term
      [1] "character"
      
      $ontology
      [1] "character"
      
      $over_represented_adj_pval
      [1] "numeric"
      
      $under_represented_adj_pval
      [1] "numeric"
      
      $term_short
      [1] "character"
      

# estimate_go_overrep() works

    Code
      lapply(out, class)
    Output
      $category
      [1] "character"
      
      $over_represented_pvalue
      [1] "numeric"
      
      $under_represented_pvalue
      [1] "numeric"
      
      $numDEInCat
      [1] "integer"
      
      $numInCat
      [1] "integer"
      
      $term
      [1] "character"
      
      $ontology
      [1] "character"
      
      $over_represented_adj_pval
      [1] "numeric"
      
      $under_represented_adj_pval
      [1] "numeric"
      
      $term_short
      [1] "character"
      
      $adj_overrep
      [1] "numeric"
      

---

    Code
      lapply(out2, class)
    Output
      $category
      [1] "character"
      
      $over_represented_pvalue
      [1] "numeric"
      
      $under_represented_pvalue
      [1] "numeric"
      
      $numDEInCat
      [1] "integer"
      
      $numInCat
      [1] "integer"
      
      $term
      [1] "character"
      
      $ontology
      [1] "character"
      
      $over_represented_adj_pval
      [1] "numeric"
      
      $under_represented_adj_pval
      [1] "numeric"
      
      $term_short
      [1] "character"
      
      $adj_overrep
      [1] "numeric"
      

# remove_redundant_go() works

    Code
      lapply(out, class)
    Output
      $category
      [1] "character"
      
      $over_represented_pvalue
      [1] "numeric"
      
      $under_represented_pvalue
      [1] "numeric"
      
      $numDEInCat
      [1] "integer"
      
      $numInCat
      [1] "integer"
      
      $term
      [1] "character"
      
      $ontology
      [1] "character"
      
      $over_represented_adj_pval
      [1] "numeric"
      
      $under_represented_adj_pval
      [1] "numeric"
      
      $term_short
      [1] "character"
      
      $adj_overrep
      [1] "numeric"
      

