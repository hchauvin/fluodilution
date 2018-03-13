# Generate documentation
devtools::document(roclets=c('rd', 'collate', 'namespace'))
system("rm -f inst/doc/fluodilution.pdf")
system("R CMD Rd2pdf -o inst/doc/fluodilution.pdf --no-preview ./")
devtools::build_vignettes()

# For spelling errors
sink("scripts/errormsg.txt")
print(tools::xgettext(".", verbose=TRUE))
sink(NULL)