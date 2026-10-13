###############################################################################
# EBC4008 Case 3 -- run every part in order (working directory = case3)
###############################################################################

for (f in list.files("R", pattern = "^[0-9]{2}_.*[.]R$", full.names = TRUE)) {
  message("== ", f)
  source(f, local = new.env())   # fresh env per part: parts share only results/rds/
}
