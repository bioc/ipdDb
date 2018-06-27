## Create the select() interface
##  Columns and keytypes method ##
.getDataTables <- function(con){
  tables <- RSQLite::dbListTables(con)
  tables <- tables[!tables %in% c("metadata","map_metadata","map_counts")]
  setNames(tables, tables)
}

.getLCcolnames <- function(x){
  con <- dbconn(x)
  tables <- .getDataTables(con)
  cols <- unique(unlist(lapply(tables, FUN=RSQLite::dbListFields, con=con)))
  cols <- cols[!cols %in% "_id"]
  cols
}

.cols <- function(x){
  list <- .getLCcolnames(x)
  ## dont use the primary keytype
  list <- list[!list == "GID"]
  toupper(list)
}


## Keys method ##
.deriveTableNameFromField <- function(field, x){
  con <- dbconn(x)
  tables <- .getDataTables(con)
  colTabs <- lapply(tables, FUN=RSQLite::dbListFields, con=con)
  m <- unlist2(lapply(colTabs, match, field))
  tab <- names(m)[!is.na(m)]
  if(length(tab) > 1){stop("Two fields in the source DB have the same name.")}
  if(length(tab) == 0){stop("Did not find a field in the source DB.")}
  tab
}

.keys <- function(x, keytype){
  keytype <- tolower(keytype)
  tab <- .deriveTableNameFromField(field=keytype, x)
  ## So now we know table name (tab) and field (keytype)
  sql <- paste("SELECT",keytype,"FROM",tab)
  res <- RSQLite::dbGetQuery(dbconn(x), sql)
  as.character(res[!is.na(res)])
}


.testForValidKeys <- function(x, keys, keytype){
  if (!is.character(keys)){
    stop("'keys' must be a character vector")
  }
  if (length(keys) == 0L) {
    return()
  }
  ktKeys <- keys(x, keytype)

  if(!(any(ktKeys %in% keys))){
    msg <- paste0("None of the keys entered are valid keys for '",keytype,
                  "'. Use the keys method to see a listing of valid arguments.")
    stop(msg) ## later when things are better, demote this to a warning()
    return(FALSE)
  }
  return(TRUE)
}

## select methods
##
.select <- function(x, keys=NULL, columns=NULL, keytype){
  keytype <- tolower(keytype)
  columns <- tolower(columns)
  .testForValidKeys(x, keys, keytype)
  ## 1st pool all the fields we need to extract
  fields <- unique(c(columns, keytype))
  ## Then get the tables to go with each one.
  tabs <- vapply(fields, .deriveTableNameFromField, x=x,
                 FUN.VALUE = character(1))
  ## make fully qualified fields of these tabs (the ones we want to extract)
  f.fields <- paste(tabs, fields, sep=".")

  ## Make non-redundant list of tables to visit
  nrTabs <- unique(tabs)
  ## Now join to each table
  for(i in seq_along(nrTabs)){
    if(i==1){
      sql <- paste("SELECT ",paste(f.fields, collapse=","),
                   " FROM",tabs[1])
    }else{
      ## IF we see c.probes in nrTabs[i], it means we have to
      ## use gene_id instead.
      if("c.probes" %in% nrTabs[i]){
        sql <- c(sql, paste("LEFT JOIN ",nrTabs[i],"USING (GID)"))
      }else{
        sql <- c(sql, paste("LEFT JOIN ",nrTabs[i],"USING (_id)"))
      }
    }
  }
  sql <- paste(sql, collapse=" ")
  ## add the where clause
  fullKeytype <- tabs[names(tabs)==keytype]
  fullKeytype <- paste(fullKeytype, names(fullKeytype), sep=".")
  strKeys <- paste0('"',keys,'"',collapse=",")
  where <- paste("WHERE ",fullKeytype,"IN (",strKeys,")" )
  sql <- paste(sql, where)
  ## then call that
  res <- RSQLite::dbGetQuery(dbconn(x), sql)
  ## cleanup and re-organize
  resort_base(res, keys, jointype=keytype, fields)
}