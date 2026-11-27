# Code release header: computes, for each mCRPC unit, the fraction of its cells with AR or KLK3 detected together with CHGA, CHGB or SYP (co-detection). Methods: Luminal subgroup annotation.
# Input: an AnnData file whose X holds raw UMI counts in CSR form, with genes indexed by symbol (such as the deposited luminal_cells_integrated.h5ad), and a cell table with columns cell_id and unit_id (not included).
# Output: ar_ne_codetection.csv, which R_rules_v9.R reads: per unit, n_cells, n_ar_or_klk3, n_chga_chgb_or_syp, n_both (cells with both) and codetect_frac = n_both / n_cells.
# Detected means 1 or more UMIs. A cell is counted if its obs disease_state is mCRPC and the cell table gives it a unit.
.libPaths(c("<PROJECT_DIR>/redo_env_v1/Rlib_R433_overlay","<R_LIBRARY>",.Library))
suppressPackageStartupMessages({library(data.table);library(rhdf5)})
H5 <- "<input: AnnData file with raw counts in X>"; CT <- "<input: cell table with columns cell_id and unit_id (the units of the labeling rule)>"; OUT <- "<output folder>"
AR_SIDE <- c("AR","KLK3"); NE_SIDE <- c("CHGA","CHGB","SYP"); CHUNK <- 5000L
rd <- function(p,...) as.vector(h5read(H5,p,bit64conversion="double",...))
stopifnot(h5readAttributes(H5,"X")[["encoding-type"]]=="csr_matrix")
cells <- rd("obs/_index"); genes <- rd("var/_index"); ds <- rd("obs/disease_state/categories")[rd("obs/disease_state/codes")+1L]
ga <- match(AR_SIDE,genes)-1L; gn <- match(NE_SIDE,genes)-1L; stopifnot(!anyNA(c(ga,gn)))
CTb <- fread(CT,select=c("cell_id","unit_id")); stopifnot(!anyDuplicated(CTb$cell_id))
ip <- rd("X/indptr"); rows <- which(ds=="mCRPC" & cells %in% CTb$cell_id); has_a <- logical(length(cells)); has_n <- has_a
## runs of consecutive rows, at most CHUNK rows per read; ix = 0-based gene index, rw = row of each stored count
for (r in split(rows,cumsum(c(1L,diff(rows)!=1L)))) for (s in split(r,ceiling(seq_along(r)/CHUNK))) {
  e <- s[length(s)]; n <- ip[e+1L]-ip[s[1]]; if (n==0) next
  ix <- rd("X/indices",start=ip[s[1]]+1,count=n); x <- rd("X/data",start=ip[s[1]]+1,count=n); rw <- rep(s,diff(ip[c(s,e+1L)]))
  has_a[unique(rw[x>=1 & ix %in% ga])] <- TRUE; has_n[unique(rw[x>=1 & ix %in% gn])] <- TRUE }
D <- data.table(unit=CTb$unit_id[match(cells[rows],CTb$cell_id)],a=has_a[rows],n=has_n[rows])
O <- D[,.(n_cells=.N,n_ar_or_klk3=sum(a),n_chga_chgb_or_syp=sum(n),n_both=sum(a & n)),by=unit][order(unit)]; O[,codetect_frac:=n_both/n_cells]
cat("mCRPC cells",length(rows),"units",nrow(O),"\n")
f <- file.path(OUT,"ar_ne_codetection.csv"); tmp <- paste0(f,".tmp"); fwrite(O,tmp); stopifnot(file.size(tmp)>0,nrow(fread(tmp))==nrow(O),file.rename(tmp,f))
