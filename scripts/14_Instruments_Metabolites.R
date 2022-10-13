diff_metabolites <- readr::read_tsv("../data/diff_metabolites.tsv")

ids_to_query <- unlist(diff_metabolites, use.names = FALSE)

metabolite_ins <- ieugwasr::tophits(ids_to_query, clump = 0)

readr::write_tsv(metabolite_ins, "../data/metabolite_ins.tsv")
