# Local Shiny entrypoint. Open the project directory with shiny::runApp(".").
for (path in list.files("R",pattern="\\.R$",full.names=TRUE)) source(path,encoding="UTF-8")
pm_wizard_app(template=normalizePath(file.path("inst","templates","report_template.html"),winslash="/",mustWork=TRUE))
