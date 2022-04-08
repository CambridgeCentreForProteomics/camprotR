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

# remove_redundant_go() works

    Code
      go_res %>% filter(grepl("immune", term)) %>% remove_redundant_go() %>% select(
        category, term, adj_overrep) %>% mutate(adj_overrep = round(adj_overrep,
        digits = 3))
    Message <simpleMessage>
      'select()' returned 1:1 mapping between keys and columns
    Output
          category
      1 GO:0002220
      2 GO:0002252
                                                                             term
      1 innate immune response activating cell surface receptor signaling pathway
      2                                                   immune effector process
        adj_overrep
      1      12.343
      2       2.081

