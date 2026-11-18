#For updates and bug fixes, please see the actively maintained version of Pb210Modeler:
#https://github.com/ram-24-mp/Pb210Modeler

# package start up
script_version <- "1.2"

required_packages <- c(
  "jsonlite", "readxl", "writexl", "ggplot2",
  "scales", "segmented", "dplyr"
)

min_versions <- c(
  rlang = "1.1.7",
  vctrs = "0.6.5"
)

get_latest_r_version <- function() {
  urls <- c(
    "https://cran.r-project.org/bin/windows/base/release.htm",
    "https://cran.r-project.org/bin/macosx/",
    "https://cran.r-project.org/"
  )
  
  for (u in urls) {
    txt <- tryCatch(readLines(u, warn = FALSE), error = function(e) NULL)
    if (is.null(txt)) next
    x <- paste(txt, collapse = "\n")
    
    m <- regmatches(
      x,
      gregexpr("\\bR[- ]([0-9]+\\.[0-9]+\\.[0-9]+)\\b", x, perl = TRUE)
    )[[1]]
    
    if (length(m) == 0) next
    
    vers <- sub("^R[- ]", "", m)
    vers <- unique(vers[nzchar(vers)])
    if (length(vers) == 0) next
    
    return(as.character(sort(numeric_version(vers), decreasing = TRUE)[1]))
  }
  
  NA_character_
}

check_r_version_current <- function() {
  current_r <- as.character(getRversion())
  latest_r <- get_latest_r_version()
  
  if (is.na(latest_r) || !nzchar(latest_r)) {
    message("Could not check whether R is up to date.")
    return(invisible(FALSE))
  }
  
  if (utils::compareVersion(current_r, latest_r) < 0) {
    stop(
      paste0(
        "Your R version is out of date (current: ", current_r,
        ", latest: ", latest_r, ").\n",
        "Please update R with install.packages('installr'); installr::updateR(), restart R, and run the script again."
      ),
      call. = FALSE
    )
  }
  
  message("R is up to date.")
  invisible(TRUE)
}

check_installed_version <- function(pkg, min_version) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(FALSE)
  utils::compareVersion(as.character(utils::packageVersion(pkg)), min_version) >= 0
}

install_missing_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  
  if (length(missing) == 0) return(invisible(FALSE))
  
  install.packages(missing, dependencies = TRUE)
  
  stop(
    paste0(
      "Installed missing packages: ",
      paste(missing, collapse = ", "),
      "\nPlease restart R and run the script again."
    ),
    call. = FALSE
  )
}

check_outdated_packages <- function(version_requirements) {
  outdated <- names(version_requirements)[
    !vapply(
      names(version_requirements),
      function(pkg) check_installed_version(pkg, version_requirements[[pkg]]),
      logical(1)
    )
  ]
  
  if (length(outdated) > 0) {
    stop(
      paste0(
        "These packages are too old and must be updated in a fresh R session: ",
        paste(outdated, collapse = ", "),
        "\nRun:\n",
        "install.packages(c(",
        paste(sprintf('"%s"', outdated), collapse = ", "),
        "), dependencies = TRUE)\n",
        "Then restart R and run the script again."
      ),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}

load_required_packages <- function(pkgs) {
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  invisible(TRUE)
}

check_for_script_update <- function(current_version) {
  latest_version <- tryCatch({
    x <- jsonlite::fromJSON(
      "https://api.github.com/repos/ram-24-mp/Pb210Modeler/releases/latest"
    )
    sub("^v", "", trimws(x$tag_name))
  }, error = function(e) NA_character_)
  
  if (is.na(latest_version) || latest_version == "") {
    message("Could not check for script updates.")
    return(invisible(FALSE))
  }
  
  if (utils::compareVersion(latest_version, current_version) > 0) {
    message(sprintf(
      "A newer version is available: %s (current: %s)",
      latest_version, current_version
    ))
    return(invisible(TRUE))
  }
  
  message("Script is up to date.")
  invisible(FALSE)
}

check_r_version_current()

if ("rlang" %in% loadedNamespaces()) {
  loaded_rlang <- as.character(getNamespaceVersion("rlang"))
  if (utils::compareVersion(loaded_rlang, min_versions[["rlang"]]) < 0) {
    stop(
      paste0(
        "An outdated 'rlang' namespace is already loaded in this R session (",
        loaded_rlang,
        ").\n",
        "Restart R before running this script.\n",
        "Do not source this script again in the same session after updating packages.\n",
        "Run install.packages(c(\"rlang\", \"vctrs\", \"ggplot2\", \"dplyr\"), dependencies = TRUE)"
      ),
      call. = FALSE
    )
  }
}

install_missing_packages(required_packages)
check_outdated_packages(min_versions)
load_required_packages(required_packages)
check_for_script_update(script_version)
# clean workspace
rm(list = ls(all.names = TRUE))

dynamic_scale <- function(table, column_name, n_breaks = 10,
                          axis = "x", reverse = FALSE) {
  
  col_data <- table[[column_name]]
  
  # Remove NA / non-finite values
  col_data <- col_data[is.finite(col_data)]
  
  # Handle empty or invalid data
  if (length(col_data) == 0) {
    stop("No finite values found in column: ", column_name)
  }
  
  # Determine the range
  data_range <- range(col_data)
  
  # Generate "nice" breaks
  brks <- scales::breaks_pretty(n = n_breaks)(data_range)
  
  # Remove duplicate/non-increasing breaks
  brks <- brks[is.finite(brks)]
  brks <- unique(brks)
  
  if (length(brks) < 2) {
    brks <- data_range
  }
  
  # Determine spacing between breaks
  diffs <- diff(brks)
  diffs <- diffs[diffs > 0]
  
  if (length(diffs) == 0) {
    accuracy <- 1
  } else {
    step <- min(diffs)
    accuracy <- 10^floor(log10(step))
  }
  
  # Return base-R breaks if requested
  if (axis == "base") {
    if (reverse) {
      brks <- rev(brks)
    }
    
    return(list(
      breaks = brks,
      accuracy = accuracy
    ))
  }
  
  # Return ggplot scale
  if (axis == "x") {
    return(
      ggplot2::scale_x_continuous(
        breaks = brks,
        labels = scales::label_number(accuracy = accuracy),
        trans = if (reverse) "reverse" else "identity"
      )
    )
  }
  
  if (axis == "y") {
    return(
      ggplot2::scale_y_continuous(
        breaks = brks,
        labels = scales::label_number(accuracy = accuracy),
        trans = if (reverse) "reverse" else "identity"
      )
    )
  }
  
  stop("axis must be 'x', 'y', or 'base'")
}

decimal_year <- function(date) {
  
  year_start <- as.Date(paste0(format(date, "%Y"), "-01-01"))
  
  next_year_start <- as.Date(paste0(
    as.integer(format(date, "%Y")) + 1, "-01-01"
  ))
  
  as.numeric(format(date, "%Y")) +
    as.numeric(date - year_start) /
    as.numeric(next_year_start - year_start)
}

# make a numeric date
day=as.numeric(readline(prompt = "Enter Day "))
month=as.numeric(readline(prompt = "Enter Month "))
year=as.numeric(readline(prompt = "Enter Year "))

# Construct actual calendar date
collection_date = as.Date(
  sprintf("%04d-%02d-%02d", year, month, day)
)
numeric_date=decimal_year(collection_date)

interval_constant=readline(prompt = "Is your interval thickness constant (True or False) ")
# core dimensions, in cm
if (interval_constant==TRUE){
  print ("Diameter 10.0 +/- 0.1 cm, Interval Thickness 1 cm")
dim_defaults=readline(prompt = "Accept Core Dimension Defaults (True or False) ")
if (dim_defaults==TRUE){
  internal_diameter=10.0
  internal_diameter_uncer=0.1
  surface_area=pi*( internal_diameter/2)^2
  surface_area_uncer=pi*internal_diameter*internal_diameter_uncer/2
  interval_thickness=1
} else {
  internal_diameter=as.numeric(readline(prompt = "Enter Diameter (default 10.0) "))
  internal_diameter_uncer=as.numeric(readline(prompt = "Enter Uncertainty (default 0.1) "))
  surface_area=pi*( internal_diameter/2)^2
  surface_area_uncer=pi*internal_diameter*internal_diameter_uncer/2
  interval_thickness=as.numeric(readline(prompt = "Enter Interval Thickness (default 1) "))
}
core_length=as.numeric(readline(prompt = "Enter Core Length "))

# Pb-210 constants: Source DDEP, 2010
half_life= 22.23
half_life_uncer=0.12
decay_const=log(2)/half_life
decay_const_uncer= decay_const*(0.12/22.23*100)/100
# DDEP year convention
days_per_year <- 365.242198

# masses
# dry bulk density g cm^-3 & corresponding uncertainty columns are provided by user. Data must be entered with empty row between each entry
mass_table=array(data=NA, dim=c((((core_length*2)/interval_thickness)+1), 11), dimnames=list(NULL, c("Top of Interval z(i) cm", "Mid Depth zi cm","Interval Thickness delta zi cm","Dry Bulk Density DBD gcm^-3","Dry Bulk Density Uncertainty u(DBD) gcm^-3","Aerial Dry Mass delta mi/S gcm^-2","Aerial Dry Mass Uncertainty u(delta mi/S) gcm^-2","Cumulative Mass Depth Layer m(i) gcm^-2","Cumulative Mass Depth Layer Uncertainty u(m(i)) gcm^-2","Mass Depth Section mi gcm^-2","Mass Depth Section Uncertainty u(mi) gcm^-2")))
dimnames(mass_table)[[2]] <- gsub("delta", "\u0394", dimnames(mass_table)[[2]])
#populate z information
for (i in 1:nrow(mass_table)){
  if (i==1){
    mass_table[i,1]=0
  } else { 
    if (i %% 2 == 1) {
      mass_table[i,1]=mass_table[i-2,1]+interval_thickness
    } else {
      mass_table[i,2]=mass_table[i-1,1]+interval_thickness/2
      mass_table[i,3]=interval_thickness
    }
  }
}
}else{ #if interval thickness is variable
  how_many_thick=as.numeric(readline(prompt = "How many different interval thicknesses are used in the core? "))
  interval_transitions_list=c()
  for (i in 1:(how_many_thick-1)){
      interval_transitions_list[i]=as.numeric(readline(prompt = paste("Enter interval thickness transition depth ", i, ": ")))
  }
  interval_thicknesses_list=c()
  for (i in 1:how_many_thick){
    interval_thicknesses_list[i]=as.numeric(readline(prompt = paste("Enter interval thickness ", i, ": ")))
  }
  print ("Diameter 10.0 +/- 0.1 cm")
  dim_defaults=readline(prompt = "Accept Core Dimension Defaults (True or False) ")
  if (dim_defaults==TRUE){
    internal_diameter=10.0
    internal_diameter_uncer=0.1
    surface_area=pi*( internal_diameter/2)^2
    surface_area_uncer=pi*internal_diameter*internal_diameter_uncer/2
    
  } else {
    internal_diameter=as.numeric(readline(prompt = "Enter Diameter (default 10.0) "))
    internal_diameter_uncer=as.numeric(readline(prompt = "Enter Uncertainty (default 0.1) "))
    surface_area=pi*( internal_diameter/2)^2
    surface_area_uncer=pi*internal_diameter*internal_diameter_uncer/2
    
  }
  interval_thickness=interval_thicknesses_list[1]
  core_length=as.numeric(readline(prompt = "Enter Core Length "))
  
  # Pb-210 constants: Source DDEP, 2010
  half_life= 22.23
  half_life_uncer=0.12
  decay_const=log(2)/half_life
  decay_const_uncer= decay_const*(0.12/22.23*100)/100
  # DDEP year convention
  days_per_year <- 365.242198
  
  
  
  interval_thicknesses <- unlist(interval_thicknesses_list)
  interval_transitions <- unlist(interval_transitions_list)
  
  # Segment boundaries: from 0, through transitions, to core_length
  bounds <- c(0, interval_transitions, core_length)
  
  # Length of each depth segment
  segment_lengths <- diff(bounds)
  
  # Number of intervals in each segment (rounded up)
  intervals_per_segment <- ceiling(segment_lengths / interval_thicknesses)
  
  # Base number of intervals
  n_intervals <- sum(intervals_per_segment)
  
  # Total rows: one for mid-depth, one for boundaries, aligned in depth order
  mass_table_rows <- 2 * n_intervals+1
  
  mass_table=array(data=NA, dim=c(mass_table_rows, 11), dimnames=list(NULL, c("Top of Interval z(i) cm", "Mid Depth zi cm","Interval Thickness delta zi cm","Dry Bulk Density DBD gcm^-3","Dry Bulk Density Uncertainty u(DBD) gcm^-3","Aerial Dry Mass delta mi/S gcm^-2","Aerial Dry Mass Uncertainty u(delta mi/S) gcm^-2","Cumulative Mass Depth Layer m(i) gcm^-2","Cumulative Mass Depth Layer Uncertainty u(m(i)) gcm^-2","Mass Depth Section mi gcm^-2","Mass Depth Section Uncertainty u(mi) gcm^-2")))
  colnames(mass_table) <- gsub("delta", "\u0394", colnames(mass_table))
  # 1) Build boundary depths
  boundary_depths <- numeric(0)
  
  for (i in seq_along(interval_thicknesses)) {
    t     <- interval_thicknesses[i]
    start <- bounds[i]
    end   <- bounds[i + 1]
    
    seg_bounds <- seq(from = start, to = end - 1e-9, by = t)
    if (tail(seg_bounds, 1) < end) {
      seg_bounds <- c(seg_bounds, end)
    }
    
    if (length(boundary_depths) > 0) {
      seg_bounds <- seg_bounds[-1]
    }
    
    boundary_depths <- c(boundary_depths, seg_bounds)
  }
  
  stopifnot(length(boundary_depths) == n_intervals + 1)
  
  # 2) Populate ONLY the first column of the existing mass_table
  #    (assumes mass_table already has mass_table_rows rows)
  boundary_rows <- seq(1, by = 2, length.out = length(boundary_depths))
  mass_table[boundary_rows, 1] <- boundary_depths
  
  # rest of z data
  for (i in seq(2, nrow(mass_table), by = 2)) {  
    mass_table[i,2]=(mass_table[(i-1),1]+mass_table[(i+1),1])/2
    mass_table[i,3]=mass_table[(i+1),1]-mass_table[(i-1),1]
  }
}

#load in mass information
# Define the path to your Excel file
print("select mass data")
file_path=file.choose()
# Read the Excel file
mass_data <- read_excel(file_path, col_names = TRUE, col_types = "numeric", .name_repair = "unique_quiet")
mass_data=as.data.frame(mass_data)

for (i in 1:nrow(mass_table)){
  if (i %% 2 == 0){
    mass_table[i,4]=mass_data[(i/2),1]
    mass_table[i,5]=mass_data[(i/2),2]
  }
}

original_section_DBD <- !is.na(
  mass_table[, "Dry Bulk Density DBD gcm^-3"]
)

original_section_DBD_evens=rep(NA,10)
for (i in 1:10){
  original_section_DBD_evens[i]=original_section_DBD[i*2]
}

# Preserve original measured values and uncertainties
original_DBD <- mass_table[, "Dry Bulk Density DBD gcm^-3"]
original_uDBD <- mass_table[, "Dry Bulk Density Uncertainty u(DBD) gcm^-3"]

# Now interpolate any missing values
for (i in 1:nrow(mass_table)){
  if (i %% 2 == 0){
    if (is.na(original_DBD[i])){
      
      # Find the previous section with an original concentration value
      j = i - 2
      while (j >= 2 && is.na(original_DBD[j])){
        j = j - 2
      }
      
      # Find the next section with an original DBD value
      k = i + 2
      while (k <= nrow(mass_table) && is.na(original_DBD[k])){
        k = k + 2
      }
      
      # Interpolate only if there is an original valid section on both sides
      if (j >= 2 && k <= nrow(mass_table)){
        
        fraction =
          (mass_table[i,"Mid Depth zi cm"] -
             mass_table[j,"Mid Depth zi cm"]) /
          (mass_table[k,"Mid Depth zi cm"] -
             mass_table[j,"Mid Depth zi cm"])
        
        mass_table[i,"Dry Bulk Density DBD gcm^-3"] =
          original_DBD[j] +
          (original_DBD[k] - original_DBD[j]) * fraction
        
        mass_table[i,"Dry Bulk Density Uncertainty u(DBD) gcm^-3"] =
          sqrt((1-fraction)^2 * original_uDBD[j]^2 +
                 fraction^2 * original_uDBD[k]^2)
      }
    }
  }
}
# aerial mass calculations
mass_table[,6]=mass_table[,3]*mass_table[,4]
mass_table[,7]=mass_table[,5]*mass_table[,3]
# cumulative mass depth
for (i in 1:nrow(mass_table)){
  if (i==1){
    mass_table[i,8]=0
  } else {
    if (i %% 2 == 1) {
      mass_table[i,8]=mass_table[i-2,8]+ mass_table[i-1,6]
    }
  }
}
# cumulative mass depth uncertainty
for (i in 1:nrow(mass_table)){
  if (i==1){
    mass_table[i,9]= 0
  } else{
    if (i %% 2 == 1) {
      mass_table[i,9]= sqrt(sum(mass_table[(i-1),7]^2, mass_table[(i-2),9]^2))
    }
  }
}
# mass depth section 
for (i in 1:nrow(mass_table)-1){
  if (i %% 2 == 0) {
    mass_table[i,10]=mass_table[i,6]*0.5+mass_table[(i-1),8]
  }
}
# mass depth section uncertainty 
for (i in 1:nrow(mass_table)-1){
  if (i %% 2 == 0) {
    mass_table[i,11]= sqrt(mass_table[(i-1),9]^2+(0.5*mass_table[i,7])^2)
  }
}
# create the concentrations table and begin preliminary population
concentrations_table=array(data=NA, dim=c((nrow(mass_table)-1), 15), dimnames=list(NULL, c("Mid Depth zi cm","Total Pb-210 dpm/g","Total Pb-210 Uncertainty u(Pb-210) dpm/g","Supported Pb-210 dpm/g","Supported Pb-210 Uncertainty u(sPb-210) dpm/g","Excess Pb-210 Ci dpm/g","Excess Pb-210 Uncertainty u(Ci) dpm/g","Excess Pb-210 C(i) dpm/g","Excess Pb-210 Uncertainty u(C(i)) dpm/g","Inventory delta Ai dpm/cm^2","Inventory Uncertainty u(delta Ai) dpm/cm^2","Mass Depth mi g/cm^2","ln(Ci)","u(ln(Ci))","Cumulative Mass Depth Layer m(i) gcm^-2")))
colnames(concentrations_table) <- gsub("delta", "\u0394", colnames(concentrations_table))
concentrations_table[,1]=mass_table[1:nrow(concentrations_table),2]
concentrations_table[,12]=mass_table[1:nrow(concentrations_table),10]
concentrations_table[,15]=mass_table[1:nrow(concentrations_table),8]
# load in Pb-210 data, assign Alpha or Gamma methods, and fill activity data
alpha_or_gamma=readline(prompt = "Do you wish to enter Alpha data (True or False) ")
if (alpha_or_gamma==TRUE){
  print("select Alpha data")
  file_path2=file.choose()
  pb210_data <- read_excel(file_path2, col_names = TRUE, col_types = "numeric", .name_repair = "unique_quiet")
  pb210_data=as.data.frame(pb210_data)
  for (i in 1:nrow(concentrations_table)){
    if (i %% 2 == 0){
      concentrations_table[i,2]=pb210_data[(i/2),1]
      concentrations_table[i,3]=pb210_data[(i/2),2]
    }
  }
  
  
  accept_background_determination=FALSE
  
  #start of background auto detection
  background_table=matrix(nrow=nrow(pb210_data), ncol=3)
  background_table=as.data.frame(background_table)
  pb210_data_auto_bkrnd=array(data=NA, dim=c(nrow(mass_table)/2,3), dimnames=list(NULL, c("y","x", "x_uncertainty")))
  for (i in 1:nrow(concentrations_table)){
    if (i%%2==0){
      background_table[i/2,1]=concentrations_table[i,1]
      background_table[i/2,2]=concentrations_table[i,2]
      background_table[i/2,3]=concentrations_table[i,3]
    }
  }
  
  pb210_data_auto_bkrng=pb210_data
  pb210_data_auto_bkrng=as.data.frame(pb210_data_auto_bkrng)
  pb210_data_auto_bkrng$X.2=na.omit(concentrations_table[,1])
  pb210_data_auto_bkrng=as.data.frame(pb210_data_auto_bkrng)
  
  #perform analysis
  detect_pb210_background <- function(x, y,
                                      hmin = 3,
                                      slope_tol = NULL,      # dpm/g per cm; if NULL, data-driven
                                      alpha_flat = 0.10,     # bottom slope == 0 (want large p)
                                      alpha_diff = 0.05,     # slope differs above vs below (want small p)
                                      min_drop = NULL,       # optional: require mean(top)-mean(bottom) >= min_drop
                                      prefer_deep = TRUE) {
    
    x <- as.numeric(x); y <- as.numeric(y)
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]; y <- y[ok]
    
    o <- order(x)
    x <- x[o]; y <- y[o]
    
    n <- length(y)
    if (n < 2 * hmin + 1) {
      return(list(present = FALSE, idx = integer(0), reason = "too_few_points"))
    }
    
    # Default slope tolerance:
    # take 10% of the absolute global slope, but never smaller than a noise-floor-based value.
    if (is.null(slope_tol)) {
      b_global <- coef(lm(y ~ x))[2]
      # noise-floor: typical activity change over 1 cm due to scatter
      noise_floor <- 0.25 * stats::mad(diff(y), constant = 1)  # dpm/g per step (~cm)
      slope_tol <- max(0.1 * abs(b_global), noise_floor)
      if (!is.finite(slope_tol) || slope_tol <= 0) slope_tol <- noise_floor
    }
    
    # Optional mean-drop requirement default: 1 MAD of y
    if (is.null(min_drop)) {
      min_drop <- stats::mad(y, constant = 1)
      if (!is.finite(min_drop) || min_drop < 0) min_drop <- 0
    }
    
    j_grid <- (hmin + 1):(n - hmin)
    
    rows <- lapply(j_grid, function(j) {
      top <- 1:(j - 1)
      bot <- j:n
      
      fit_top <- lm(y[top] ~ x[top])
      fit_bot <- lm(y[bot] ~ x[bot])
      
      b_top <- unname(coef(fit_top)[2])
      b_bot <- unname(coef(fit_bot)[2])
      
      se_top <- summary(fit_top)$coefficients[2, 2]
      se_bot <- summary(fit_bot)$coefficients[2, 2]
      
      df_top <- length(top) - 2
      df_bot <- length(bot) - 2
      
      # bottom slope == 0 test (want NOT significant)
      t0 <- b_bot / se_bot
      p_flat <- 2 * pt(abs(t0), df = df_bot, lower.tail = FALSE)
      
      # slope difference test (want significant)
      tdiff <- (b_top - b_bot) / sqrt(se_top^2 + se_bot^2)
      df_diff <- (se_top^2 + se_bot^2)^2 / ((se_top^4 / df_top) + (se_bot^4 / df_bot))
      p_diff <- 2 * pt(abs(tdiff), df = df_diff, lower.tail = FALSE)
      
      mu_top <- mean(y[top])
      mu_bot <- mean(y[bot])
      drop <- mu_top - mu_bot
      
      data.frame(
        j = j,
        x_at_start = x[j],
        n_top = length(top),
        n_bot = length(bot),
        slope_top = b_top,
        slope_bot = b_bot,
        p_flat = p_flat,
        p_diff = p_diff,
        drop = drop,
        abs_slope_bot = abs(b_bot)
      )
    })
    
    tab <- do.call(rbind, rows)
    
    keep <- with(tab,
                 abs_slope_bot <= slope_tol &
                   p_diff <= alpha_diff &
                   drop >= min_drop
    )
    
    if (!any(keep)) {
      return(list(
        present = FALSE,
        idx = integer(0),
        slope_tol = slope_tol,
        min_drop = min_drop,
        candidates = tab
      ))
    }
    
    cand <- tab[keep, , drop = FALSE]
    best <- cand[order(cand$j), ][1, , drop = FALSE]   # earliest passing j
    
    list(
      present = TRUE,
      idx = best$j,            # index in x/y after sorting by depth
      x_at_start = best$x_at_start,
      slope_bot = best$slope_bot,
      slope_tol = slope_tol,
      min_drop = min_drop,
      best_row = best,
      candidates = tab
    )
  }
  
  # Use with your columns (y = dpm/g, x = cm)
  y <- pb210_data_auto_bkrng[, 1]
  x <- pb210_data_auto_bkrng[, 3]
  res <- detect_pb210_background(x, y, hmin = 3)
  res$present
  res$idx
  res$x_at_start
  res$slope_bot
  if (res$present==TRUE) {
    suggested_background_index <- 2*res$idx
    cat("Suggested background activity begins at index:", suggested_background_index, "\n")
    print("Arithmetic Mean & Standard Error")
    print("Pb-210 Background Activity (dpm/g)")
    print(mean(pb210_data_auto_bkrng[res$idx:length(pb210_data_auto_bkrng[,1]), 1]))
    print("Pb-210 Background Activity Uncertainty (dpm/g)")
    length_finder_1=length(pb210_data_auto_bkrng[res$idx:length(pb210_data_auto_bkrng[,1]), 1])
    print(sd(pb210_data_auto_bkrng[res$idx:length(pb210_data_auto_bkrng[,1]), 1])/sqrt(length_finder_1))
  } else {
    cat("Background not detected.\n")
  }
  # end of background auto detection
  
  
  while (accept_background_determination==FALSE){
    points_to_view=as.numeric(readline(prompt = "How many datapoints do you wish to view for background determination? "))
    con_tab_rows=nrow(concentrations_table)
    points_to_view2=2*points_to_view
    start_row_view=con_tab_rows-points_to_view2
    row_name_indexs=seq(1:nrow(concentrations_table))
    row_name_indexs=row_name_indexs[start_row_view:nrow(concentrations_table)]
    row_name_indexs=list(as.character(row_name_indexs))
    background_determination_table=array(concentrations_table[start_row_view:con_tab_rows,2], dim=c(points_to_view2+1,1), dimnames=c(row_name_indexs,"Activity"))
    background_determination_table=na.omit(background_determination_table)
    attributes(background_determination_table)$na.action <- NULL
    background_determination_table_subset <- background_determination_table[2:nrow(background_determination_table), ]
    # Convert to data frame
    df <- as.data.frame(background_determination_table_subset)
    # Rename the column
    colnames(df) <- "Pb-210 Activity (dpm/g)"
    # Print the data frame
    print(df)
    background_table=array(data=NA, dim=c(nrow(mass_table)/2,3), dimnames=list(NULL, c("x","y", "x_uncertainty")))
    for (i in 1:nrow(concentrations_table)){
      if (i%%2==0){
        background_table[i/2,1]=concentrations_table[i,2]
        background_table[i/2,2]=concentrations_table[i,1]
        background_table[i/2,3]=concentrations_table[i,3]
      }
    }
    
    background_table=as.data.frame(background_table)
    background_table_plot <- ggplot(
      background_table,
      aes(x = x, y = y)
    ) +
      geom_errorbar(
        aes(
          xmin = x - x_uncertainty,
          xmax = x + x_uncertainty
        ),
        width = 0.2
      ) +
      geom_point() +
      labs(
        title = "Total Pb-210 vs Mid Depth",
        x = "Total Pb-210 (dpm/g)",
        y = "Mid Depth (cm)"
      ) +
      dynamic_scale(
        background_table,
        "x"
      ) +
      dynamic_scale(
        background_table,
        "y",
        axis = "y",
        reverse = TRUE
      ) +
      theme_minimal()
    
    print(background_table_plot)
    
    hit_background <- tolower(readline(
      prompt = "Was equilibrium depth achieved? (True or False) "
    )) == "true"
    
    if (hit_background==TRUE){
    first_background=as.numeric(readline(prompt = "At what index does background begin? "))
    background_value=mean(na.omit(concentrations_table[first_background:nrow(concentrations_table),2]))
    check_for_negs=which((concentrations_table[1:(first_background-1),2]-background_value)<0)[1]
    check_background_length=length(concentrations_table[first_background:nrow(concentrations_table)])
    test_length_finder=length(na.omit(concentrations_table[first_background:nrow(concentrations_table),2]))

    if (is.na(check_for_negs)==FALSE){
      accept_new_background=FALSE
      while(accept_new_background==FALSE){
        cat(sprintf("poor background selection detected! negative excess Pb-210 activity at index: %f\n",check_for_negs))
        accept_new_background=readline(prompt = "Accept new background determination from this index? (True) ")
        background_value=mean(na.omit(concentrations_table[check_for_negs:nrow(concentrations_table),2]))
        background_value_uncer=sd(na.omit(concentrations_table[check_for_negs:nrow(concentrations_table),2]))/sqrt(test_length_finder)
        first_background=check_for_negs
      }
    }else{
      background_value=mean(na.omit(concentrations_table[first_background:nrow(concentrations_table),2]))
      background_value_uncer=sd(na.omit(concentrations_table[first_background:nrow(concentrations_table),2]))/sqrt(test_length_finder)
    }
    if (check_background_length<5){
      print("Background selection does not contain enough points! Make a new selection")
      accept_background_determination=FALSE
    }else{
    print("Arithmetic Mean & Standard Error")
    print("Background activity (dpm/g)")
    print(background_value)
    print("Background activity uncertainty (dpm/g)")
    print(background_value_uncer)
    accept_background_determination=readline(prompt = "Accept background determination? (True or False) ")
    }
    }else{
      break
    }
  }
  if (hit_background==TRUE){
  for (i in 1:nrow(concentrations_table)){
    if (i %% 2 == 0 && is.na(concentrations_table[i,2])== FALSE) {
      concentrations_table[i,4]=background_value
      concentrations_table[i,5]=background_value_uncer
      }
    }
  
  concentrations_table=concentrations_table[1:first_background-1,]
  mass_table=mass_table[1:first_background-1,]
  }else{
    independent_background=readline(prompt = "Do you have an independent estimate of supported Pb-210 activity? (True or False) ")
    if (independent_background==TRUE){
    ind_bkrnd=as.numeric(readline(prompt = "Enter supported activity (dpm/g) "))
    ind_bkrnd_uncer=as.numeric(readline(prompt = "Enter supported activity uncertainty (dpm/g) "))
    for (i in 1:nrow(concentrations_table)){
      if (i %% 2 == 0 && is.na(concentrations_table[i,2])== FALSE) {
        concentrations_table[i,4]=ind_bkrnd
        concentrations_table[i,5]=ind_bkrnd_uncer
      }
    }
    }else{
    for (i in 1:nrow(concentrations_table)){
      if (i %% 2 == 0 && is.na(concentrations_table[i,2])== FALSE) {
          concentrations_table[i,4]=0
          concentrations_table[i,5]=0
      }
    }
    }
  }
} else {
  print("select Gamma data")
  # Define the path to your Excel file
  file_path2=file.choose()
  pb210_data <- read_excel(file_path2, col_names = TRUE, col_types = "numeric", .name_repair = "unique_quiet")
  pb210_data=as.data.frame(pb210_data)
  for (i in 1:nrow(concentrations_table)){
    if (i %% 2 == 0){
      concentrations_table[i,2]=pb210_data[(i/2),1]
      concentrations_table[i,3]=pb210_data[(i/2),2]
      concentrations_table[i,4]=pb210_data[(i/2),3]
      concentrations_table[i,5]=pb210_data[(i/2),4]
    }
  }
}


# Preserve which section Ci values were originally measured.
# These are the only section values that should receive ln(Ci).
original_section_Ci <- !is.na(concentrations_table[, "Total Pb-210 dpm/g"])

original_section_Ci_evens=rep(NA,10)
for (i in 1:10){
  original_section_Ci_evens[i]=original_section_Ci[i*2]
}
# Section ln(Ci) and u(ln(Ci)) must remain NA when
# the section Ci was not originally measured.
for (i in seq(2, nrow(concentrations_table), by = 2)) {
  if (!original_section_Ci[i]) {
    concentrations_table[i, "ln(Ci)"] <- NA_real_
    concentrations_table[i, "u(ln(Ci))"] <- NA_real_
  }
}
  self_absorb_corr=readline(prompt="Do you wish to correct total Pb-210 Activity for self absorption? (TRUE or FALSE) ")
  if (self_absorb_corr==TRUE){
    print("select attenuation coefficient data")
    self_ab_file_path <- file.choose()
    
    attenuation_data <- read_excel(
      self_ab_file_path,
      col_names = TRUE,
      .name_repair = "unique_quiet"
    )
    correction_summary <- as.data.frame(matrix(NA, nrow = nrow(attenuation_data) * 2, ncol = ncol(attenuation_data)))
    
    colnames(correction_summary) <- colnames(attenuation_data)
    
    correction_summary[seq(2, nrow(correction_summary), by = 2), ] <- attenuation_data
    correction_summary$`Original Total Activity (dpm/g)`=concentrations_table[,2]
    correction_summary$`Original Total Activity Uncertainty (dpm/g)`=concentrations_table[,3]
    correction_summary$`Corrected Total Activity (dpm/g)`= correction_summary$`Original Total Activity (dpm/g)`*(1/correction_summary[,1])
    correction_summary$`Corrected Total Activity Uncertainty (dpm/g)`=correction_summary$`Corrected Total Activity (dpm/g)` * sqrt((correction_summary[,2]/correction_summary[,1])^2+(correction_summary[,4]/correction_summary[,3])^2)
    concentrations_table[,2] <-
      correction_summary$`Corrected Total Activity (dpm/g)`
    
    concentrations_table[,3] <-
      correction_summary$`Corrected Total Activity Uncertainty (dpm/g)`
    
  }


concentrations_table[,6]=concentrations_table[,2]-concentrations_table[,4]
negative_index=which(concentrations_table[,6]<0)
if (length(negative_index) > 0) {
  concentrations_table=concentrations_table[1:(negative_index[1]-1),]
  mass_table=mass_table[1:(negative_index[1]-1),]
  original_section_Ci = original_section_Ci[1:nrow(concentrations_table)]
}
if (is.na(concentrations_table[nrow(concentrations_table),1])==TRUE){
  concentrations_table=concentrations_table[1:(nrow(concentrations_table)-1),]
  original_section_Ci = original_section_Ci[1:nrow(concentrations_table)]
}

# excess pb-210 uncertainty
for (i in 1:nrow(concentrations_table)){
  if (i %% 2 == 0) {
    concentrations_table[i,7]=sqrt(sum(concentrations_table[i,3]^2, concentrations_table[i,5]^2))
  }
}


perform_decay_corr=readline(prompt = "Do you want to decay correct excess Pb-210 activity? (TRUE or FALSE) ")
if (perform_decay_corr==TRUE){
# ==========================================
# EXCESS Pb-210 DECAY CORRECTION
# ==========================================

# Read measurement dates from Excel

print("Select Excel file containing measurement dates")

date_file_path <- file.choose()

date_data <- read_excel(
  date_file_path,
  col_names = TRUE,
  .name_repair = "unique_quiet"
)

measurement_dates_raw <- date_data[[1]]
measurement_dates_raw=measurement_dates_raw[1:(nrow(concentrations_table)/2)]

# Convert Excel dates into R Date format

if (inherits(measurement_dates_raw, "Date")) {
  
  measurement_dates <- measurement_dates_raw
  
} else if (inherits(measurement_dates_raw, "POSIXt")) {
  
  measurement_dates <- as.Date(measurement_dates_raw)
  
} else if (is.numeric(measurement_dates_raw)) {
  
  measurement_dates <- as.Date(
    measurement_dates_raw,
    origin = "1899-12-30"
  )
  
} else {
  
  measurement_dates <- as.Date(
    trimws(as.character(measurement_dates_raw)),
    format = "%m/%d/%Y"
  )
}

if (anyNA(measurement_dates)) {
  stop("One or more measurement dates could not be interpreted.")
}

# Convert dates to decimal years

measurement_numeric_year <- decimal_year(measurement_dates)

# Check number of dates against number of measured sections

measured_rows <- seq(2, nrow(concentrations_table), by = 2)

if (length(measurement_numeric_year) != length(measured_rows)) {
  stop("Number of dates does not match number of sample intervals.")
}

# Calculate elapsed time

elapsed_days = as.numeric(measurement_dates - collection_date)

elapsed_time = elapsed_days / days_per_year

if (any(elapsed_time < 0)) {
  stop("A measurement date occurs after the reference date.")
}

# Calculate decay correction factor

decay_factor <- exp(decay_const * elapsed_time)

# Preserve original excess activity and uncertainty

original_excess <- concentrations_table[measured_rows, 6]
original_excess_uncertainty <- concentrations_table[measured_rows, 7]

# Correct excess Pb-210 activity

concentrations_table[measured_rows, 6] <-
  original_excess * decay_factor

# Propagate uncertainty, including half-life uncertainty

decay_const_uncer <- decay_const *
  half_life_uncer / half_life

concentrations_table[measured_rows, 7] <- sqrt(
  (decay_factor * original_excess_uncertainty)^2 +
    (concentrations_table[measured_rows, 6] *
       elapsed_time * decay_const_uncer)^2
)

# Save date correction information

decay_summary <- data.frame(
  Measurement_Date = measurement_dates,
  Measurement_Year = measurement_numeric_year,
  Reference_Year = numeric_date,
  Elapsed_Time_Years = elapsed_time,
  Decay_Factor = decay_factor,
  Original_Excess = original_excess,
  Original_Uncertainty = original_excess_uncertainty,
  Corrected_Excess = concentrations_table[measured_rows, 6],
  Corrected_Uncertainty = concentrations_table[measured_rows, 7]
)
}



if (is.na(concentrations_table[2,6])){
# Construct the first 10 interval records
C1_table <- data.frame(
  zi = mass_table[seq(2, 20, by = 2), 2],
  mi = concentrations_table[seq(2, 20, by = 2), 12],
  ln_Ci = NA_real_
)

activity <- concentrations_table[
  seq(2, 20, by = 2),
  "Excess Pb-210 Ci dpm/g"
]

# Only take logarithms of finite, positive activities
valid_activity <- is.finite(activity) & activity > 0

C1_table$ln_Ci[valid_activity] <- log(activity[valid_activity])

# Preserve the original-measurement rule
for (j in seq_len(nrow(C1_table))) {
  if (is.na(original_section_Ci_evens[j]) ||
      !original_section_Ci_evens[j]) {
    C1_table$ln_Ci[j] <- NA_real_
  }
}

C1_model_list <- list()
model_summaries_surf_tot_activity <- list()

for (extent in 3:10) {
  
  candidate <- C1_table[seq_len(extent), ]
  
  # Identify valid observations
  valid <- is.finite(candidate[, "mi"]) &
    is.finite(candidate[, "ln_Ci"])
  
  fit_data <- candidate[valid, , drop = FALSE]
  
  # Use the same columns for every candidate model
  result <- data.frame(
    Model = paste0("Indices 1-", extent),
    N = nrow(fit_data),
    Used_indices = paste(which(valid), collapse = ", "),
    Intercept = NA_real_,
    Slope = NA_real_,
    Intercept_SE = NA_real_,
    Slope_SE = NA_real_,
    R.squared = NA_real_,
    Adj.R.squared = NA_real_,
    Status = "Insufficient valid observations",
    stringsAsFactors = FALSE
  )
  
  # Require at least 3 observations and distinct predictor values
  if (nrow(fit_data) >= 3 &&
      length(unique(fit_data$mi)) >= 2) {
    
    model <- lm(ln_Ci ~ mi, data = fit_data)
    
    C1_model_list[[as.character(extent - 2)]] <- model
    
    sm <- summary(model)
    coefs <- sm$coefficients
    
    result$Intercept <- unname(coef(model)[1])
    result$Slope <- unname(coef(model)[2])
    result$Intercept_SE <- coefs[1, "Std. Error"]
    result$Slope_SE <- coefs[2, "Std. Error"]
    result$R.squared <- sm$r.squared
    result$Adj.R.squared <- sm$adj.r.squared
    result$Status <- "Fitted"
  }
  
  model_summaries_surf_tot_activity[[as.character(extent)]] <- result
}

# Combine all candidate summaries
model_summary_table <- do.call(
  rbind,
  model_summaries_surf_tot_activity
)
model_summaries_surf_tot_activity=model_summary_table
rownames(model_summaries_surf_tot_activity) <- NULL
print("Surface activity missing! Please determine total surface activity")

activity <- concentrations_table[
  seq(2, 20, by = 2),
  "Excess Pb-210 Ci dpm/g"
]
C1_table[[3]] <- NA_real_
names(C1_table)[3] <- "Ci"
# Only take logarithms of finite, positive activities
valid_activity <- is.finite(activity) & activity > 0

C1_table$Ci[valid_activity] <- activity[valid_activity]

print(C1_table)
print(model_summaries_surf_tot_activity[,-(2:3)])
selection_complete=FALSE
while (selection_complete==FALSE){
  best_model_index=readline(prompt = "Select model (1-8) ")
  
  # Step 3: Retrieve the corresponding model
  best_model <- C1_model_list[[best_model_index]]
  m1 <- C1_table[1,2]
  
  beta <- coef(best_model)
  V <- vcov(best_model)
  
  log_C1 <- beta[1] + beta[2] * m1
  log_C1_uncer <- sqrt(
    V[1,1] +
      m1^2 * V[2,2] +
      2 * m1 * V[1,2]
  )
  
  C1 <- exp(log_C1)
  C1_uncer <- C1 * log_C1_uncer
  print("Surface Activity (dpm/g)")
  print(unname(C1))
  selection_complete=readline(prompt = "Accept model selection? (TRUE or FALSE) ")
}  
}

# Now interpolate any missing Total Pb-210 section values
# using exponential interpolation.
#Exponential interpolation:
# B = A^(1-r) * C^r
#Propagated uncertainty
#u(B) = B * sqrt[(1-r)^2 * (u(A)/A)^2 + r^2 * (u(C)/C)^2]

for (i in seq(2, nrow(concentrations_table), by = 2)){
  # surface absence check
  if (is.na(concentrations_table[2,"Excess Pb-210 Ci dpm/g"])){
    concentrations_table[2,"Excess Pb-210 Ci dpm/g"]=C1
    concentrations_table[2,"Excess Pb-210 Uncertainty u(Ci) dpm/g"]=C1_uncer
  }
  # Check whether this section is missing
  if (is.na(concentrations_table[i,"Excess Pb-210 Ci dpm/g"])){
    
    # Find the first missing section in this contiguous gap
    gap_start = i
    
    # Find the last missing section in this contiguous gap
    gap_end = i
    
    while (
      gap_end + 2 <= nrow(concentrations_table) &&
      is.na(concentrations_table[gap_end + 2,
                                 "Excess Pb-210 Ci dpm/g"])
    ){
      gap_end = gap_end + 2
    }
    
    # Find the valid section immediately above the gap
    j = gap_start - 2
    
    # Find the valid section immediately below the gap
    k = gap_end + 2
    
    # Only interpolate if there is a valid section
    # on both sides of the entire gap
    if (j >= 2 &&
        k <= nrow(concentrations_table) &&
        !is.na(concentrations_table[j,"Excess Pb-210 Ci dpm/g"]) &&
        !is.na(concentrations_table[k,"Excess Pb-210 Ci dpm/g"])){
      
      # Store the two bounding section concentrations
      A = concentrations_table[j,"Excess Pb-210 Ci dpm/g"]
      C = concentrations_table[k,"Excess Pb-210 Ci dpm/g"]
      
      # Store their uncertainties
      u_A =
        concentrations_table[j,
                             "Excess Pb-210 Uncertainty u(Ci) dpm/g"]
      
      u_C =
        concentrations_table[k,
                             "Excess Pb-210 Uncertainty u(Ci) dpm/g"]
      
      # Fill every missing section in the gap
      for (l in seq(gap_start, gap_end, by = 2)){
        
        # Fractional position of this section
        # between the two bounding valid sections
        fraction =
          (concentrations_table[l,"Mass Depth mi g/cm^2"] -
             concentrations_table[j,"Mass Depth mi g/cm^2"]) /
          (concentrations_table[k,"Mass Depth mi g/cm^2"] -
             concentrations_table[j,"Mass Depth mi g/cm^2"])
        
        concentrations_table[l,"Excess Pb-210 Ci dpm/g"] =
          A^(1-fraction) * C^fraction
        
        concentrations_table[l,
                             "Excess Pb-210 Uncertainty u(Ci) dpm/g"] =
          concentrations_table[l,"Excess Pb-210 Ci dpm/g"] *
          sqrt(
            (1-fraction)^2 * (u_A / A)^2 +
              fraction^2 * (u_C / C)^2
          )
      }
    }
  }
}


#C0CA calculations
did_measure_C0=readline(prompt = "Do you have a measured surface layer Pb-210 activity C(0)? (True or False) ")
if (did_measure_C0==TRUE){
  measured_C0=as.numeric(readline(prompt = "Enter measured C(0) dpm/g "))
  measured_C0_uncer=as.numeric(readline(prompt = "Enter measured C(0) uncertainty dpm/g "))
}else{
  C0CA_table=array(data=NA, dim=c(10,3), dimnames=list(NULL, c("zi","mi","ln(Ci)")))
  
  # Construct the first 10 interval records
  C0CA_table <- data.frame(
    zi = mass_table[seq(2, 20, by = 2), 2],
    mi = concentrations_table[seq(2, 20, by = 2), 12],
    ln_Ci = NA_real_
  )
  
  activity <- concentrations_table[
    seq(2, 20, by = 2),
    "Excess Pb-210 Ci dpm/g"
  ]
  
  # Only take logarithms of finite, positive activities
  valid_activity <- is.finite(activity) & activity > 0
  
  C0CA_table$ln_Ci[valid_activity] <- log(activity[valid_activity])
  
  # Preserve the original-measurement rule
  for (j in seq_len(nrow(C0CA_table))) {
    if (is.na(original_section_Ci_evens[j]) ||
        !original_section_Ci_evens[j]) {
      C0CA_table$ln_Ci[j] <- NA_real_
    }
  }
  
  C0CA_model_list <- list()
  model_summaries_C0CA_activity <- list()
  
  for (extent in 3:10) {
    
    candidate <- C0CA_table[seq_len(extent), ]
    
    # Identify valid observations
    valid <- is.finite(candidate[, "mi"]) &
      is.finite(candidate[, "ln_Ci"])
    
    fit_data <- candidate[valid, , drop = FALSE]
    
    # Use the same columns for every candidate model
    result <- data.frame(
      Model = paste0("Indices 1-", extent),
      N = nrow(fit_data),
      Used_indices = paste(which(valid), collapse = ", "),
      Intercept = NA_real_,
      Slope = NA_real_,
      Intercept_SE = NA_real_,
      Slope_SE = NA_real_,
      R.squared = NA_real_,
      Adj.R.squared = NA_real_,
      Status = "Insufficient valid observations",
      stringsAsFactors = FALSE
    )
    
    # Require at least 3 observations and distinct predictor values
    if (nrow(fit_data) >= 3 &&
        length(unique(fit_data$mi)) >= 2) {
      
      model <- lm(ln_Ci ~ mi, data = fit_data, na.action = na.exclude)
      
      C0CA_model_list[[as.character(extent - 2)]] <- model
      
      sm <- summary(model)
      coefs <- sm$coefficients
      
      result$Intercept <- unname(coef(model)[1])
      result$Slope <- unname(coef(model)[2])
      result$Intercept_SE <- coefs[1, "Std. Error"]
      result$Slope_SE <- coefs[2, "Std. Error"]
      result$R.squared <- sm$r.squared
      result$Adj.R.squared <- sm$adj.r.squared
      result$Status <- "Fitted"
    }
    
    model_summaries_C0CA_activity[[as.character(extent)]] <- result
  }
  
  # Combine all candidate summaries
  model_summary_table <- do.call(
    rbind,
    model_summaries_C0CA_activity
  )
  model_summaries_C0CA_activity=model_summary_table
  rownames(model_summaries_C0CA_activity) <- NULL
  print("Please determine surface layer excess Pb-210 activity")
  
  activity <- concentrations_table[
    seq(2, 20, by = 2),
    "Excess Pb-210 Ci dpm/g"
  ]
  C0CA_table[[3]] <- NA_real_
  names(C0CA_table)[3] <- "Ci"
  # Only take logarithms of finite, positive activities
  valid_activity <- is.finite(activity) & activity > 0
  
  C0CA_table$Ci[valid_activity] <- activity[valid_activity]
  
  # Preserve the original-measurement rule
  for (j in seq_len(nrow(C0CA_table))) {
    if (is.na(original_section_Ci_evens[j]) ||
        !original_section_Ci_evens[j]) {
      C0CA_table$Ci[j] <- NA_real_
    }
  }
  print(C0CA_table)
  print(model_summaries_C0CA_activity[,-(2:3)])
  selection_complete=FALSE
  while (selection_complete==FALSE){
    best_model_index=readline(prompt = "Select model (1-8) ")
    
    # Step 3: Retrieve the corresponding model
    best_model <- C0CA_model_list[[best_model_index]]
    
    C0CA <- exp(best_model$coefficients[1])
    C0CA_uncer <- C0CA*summary(best_model)$coefficients[1, 2]
    print("Surface Activity (dpm/g)")
    print(unname(C0CA))
    selection_complete=readline(prompt = "Accept model selection? (TRUE or FALSE) ")
  }  
}

# inventory, inventory uncertainty, mass depth, ln(Ci), u(ln(Ci)) calculations

for (i in seq(2, nrow(concentrations_table), by = 2)) {
  
  # Inventory can use interpolated section Ci
  concentrations_table[i,10] =
    concentrations_table[i,6] * mass_table[i,6]
  
  concentrations_table[i,11] =
    concentrations_table[i,10] *
    sqrt(
      sum(
        (concentrations_table[i,7] /
           concentrations_table[i,6])^2,
        (mass_table[i,7] /
           mass_table[i,6])^2
      )
    )
    concentrations_table[i,13] =
      log(concentrations_table[i,6])
    
    concentrations_table[i,14] =
      concentrations_table[i,7] /
      concentrations_table[i,6]
}

# C(i)
# Exponential Pb-210 concentration interpolation
# performed as linear interpolation in ln(C)

for (i in 1:nrow(concentrations_table)){
  if (i %% 2 == 1){
    
    if (i == 1){
      
      if(did_measure_C0==FALSE){
        concentrations_table[i,8]=C0CA
        concentrations_table[i,9]=C0CA_uncer
      }else{
        concentrations_table[i,8]=measured_C0
        concentrations_table[i,9]=measured_C0_uncer
      }
      
    } else {
      
      # Fractional position of layer between bounding sections
      fraction =
        (concentrations_table[i,"Cumulative Mass Depth Layer m(i) gcm^-2"] -
           concentrations_table[(i-1),"Mass Depth mi g/cm^2"]) /
        (concentrations_table[(i+1),"Mass Depth mi g/cm^2"] -
           concentrations_table[(i-1),"Mass Depth mi g/cm^2"])
      
      # Bounding section log concentrations
      ln_A = concentrations_table[(i-1),"ln(Ci)"]
      ln_C = concentrations_table[(i+1),"ln(Ci)"]
      
      # Bounding section log uncertainties
      u_ln_A = concentrations_table[(i-1),"u(ln(Ci))"]
      u_ln_C = concentrations_table[(i+1),"u(ln(Ci))"]
      
      # Linear interpolation in log space
      
      concentrations_table[i,"ln(Ci)"] =
        (1 - fraction) * ln_A +
        fraction * ln_C
      
      # Propagated uncertainty in log space
      
      concentrations_table[i,"u(ln(Ci))"] =
        sqrt(
          (1 - fraction)^2 * u_ln_A^2 +
            fraction^2 * u_ln_C^2
        )
      
      # Convert log concentration back to concentration
      
      concentrations_table[i,"Excess Pb-210 C(i) dpm/g"] =
        exp(concentrations_table[i,"ln(Ci)"])
      
      # Convert log uncertainty back to concentration uncertainty
      
      concentrations_table[i,"Excess Pb-210 Uncertainty u(C(i)) dpm/g"] =
        concentrations_table[i,"Excess Pb-210 C(i) dpm/g"] *
        concentrations_table[i,"u(ln(Ci))"]
    }
  }
}
for (i in 1:nrow(concentrations_table)){
  if (original_section_Ci[i]==FALSE){
    concentrations_table[i,13]=NA
    concentrations_table[i,14]=NA
  }
}
#save off copies of background-corrected mass and concentration data, accompanying plots
# find your destination folder
print("select destination folder")
destination_folder <- rstudioapi::selectDirectory("Select Destination Folder")

# Check if a folder was selected
if (!is.null(destination_folder) && destination_folder != "") {
  # Create the new folder path
  new_folder_path <- file.path(destination_folder, "Pb210_results")
  
  # Create the new folder
  dir.create(new_folder_path, showWarnings = FALSE)
  
  # Navigate into the new directory
  setwd(new_folder_path)
  
  # Print the current working directory to confirm
  print(paste("Current working directory:", getwd()))
} else {
  print("No folder selected.")
}
#start saving process
folder_name = "Background corrected mass and concentration data"

# Create the folder in the working directory
dir.create(folder_name)

# Define the file path where you want to save the file
mass_file_path <- file.path(folder_name, "mass_table.xlsx")

# Save the table as an Excel file
write_xlsx(as.data.frame(mass_table), mass_file_path)

# Define the file path where you want to save the file
concentrations_file_path <- file.path(folder_name, "concentrations_data.xlsx")

# Save the table as an Excel file
write_xlsx(as.data.frame(concentrations_table), concentrations_file_path)

#C0CA save
write_xlsx(model_summaries_C0CA_activity, path = "Background corrected mass and concentration data/C0CA_model_summaries.xlsx")
#missing surface save

if (exists("C1_model_list")){
  write_xlsx(model_summaries_surf_tot_activity, path = "Background corrected mass and concentration data/surf_ex_activity_model_summaries.xlsx")
}

if (exists("decay_summary")){
  write_xlsx(decay_summary, path = "Background corrected mass and concentration data/decay_correction_summary.xlsx")
}

# make plots
mass_depth_table=data.frame(na.omit(mass_table[,10]), na.omit(mass_table[,11]))
colnames(mass_depth_table) <- c("Mass Depth mi (g/cm^2)", "Mass Depth Uncertainty u(mi) (g/cm^2)")
#dbd
pdf(file.path(folder_name, "dry_bulk_plot.pdf"), width = 8, height = 6)
dbd_plotting_table=array(data=NA, dim=c(nrow(mass_table)/2,3), dimnames=list(NULL, c("x","y","y_uncertainty")))
dbd_plotting_table[,1]=na.omit(mass_table[,2])
dbd_plotting_table[,2]=na.omit(mass_table[,4])
dbd_plotting_table[,3]=na.omit(mass_table[,5])
dbd_plotting_table=as.data.frame(dbd_plotting_table)
# Generate breaks dynamically

dbd_plotting_table <- as.data.frame(dbd_plotting_table)

dbd_plot <- ggplot(dbd_plotting_table, aes(x = x, y = y)) +
  geom_point() +
  geom_errorbar(
    aes(
      ymin = y - y_uncertainty,
      ymax = y + y_uncertainty
    ),
    width = 0.2
  ) +
  labs(
    title = "Dry Bulk Density vs Depth",
    x = "Mid Depth (cm)",
    y = "Dry Bulk Density (g/cm^3)"
  ) +
  dynamic_scale(
    dbd_plotting_table,
    "x",
    axis = "x"
  ) +
  dynamic_scale(
    dbd_plotting_table,
    "y",
    axis = "y"
  ) +
  theme_minimal()

print(dbd_plot)
dev.off()
print(dbd_plot)
# excess pb-210
pdf(file.path(folder_name, "excess_pb210.pdf"), width = 8, height = 6)
epb210_plotting_table=array(data=NA, dim=c(nrow(mass_table)/2,3), dimnames=list(NULL, c("x","y","x_uncertainty")))
epb210_plotting_table[,1]=na.omit(concentrations_table[,6])
epb210_plotting_table[,2]=na.omit(concentrations_table[,1])
epb210_plotting_table[,3]=na.omit(concentrations_table[,7])
epb210_plotting_table=as.data.frame(epb210_plotting_table)
epb210_plot <- ggplot(
  epb210_plotting_table,
  aes(x = x, y = y)
) +
  geom_point() +
  geom_errorbar(
    aes(
      xmin = x - x_uncertainty,
      xmax = x + x_uncertainty
    ),
    width = 0.2
  ) +
  labs(
    title = "Excess Pb-210 vs Depth",
    x = "Excess Pb-210 (dpm/g)",
    y = "Mid Depth (cm)"
  ) +
  dynamic_scale(
    epb210_plotting_table,
    "x"
  ) +
  dynamic_scale(
    epb210_plotting_table,
    "y",
    axis = "y",
    reverse = TRUE
  ) +
  theme_minimal()

print(epb210_plot)
dev.off()
# supported pb-210
pdf(file.path(folder_name, "supported_pb210.pdf"), width = 8, height = 6)
spb210_plotting_table=array(data=NA, dim=c(nrow(mass_table)/2,3), dimnames=list(NULL, c("x","y","x_uncertainty")))
spb210_plotting_table[,1]=concentrations_table[seq(2, nrow(concentrations_table), by=2),4]
spb210_plotting_table[,2]=concentrations_table[seq(2, nrow(concentrations_table), by=2),1]
spb210_plotting_table[,3]=concentrations_table[seq(2, nrow(concentrations_table), by=2),5]
spb210_plotting_table=as.data.frame(spb210_plotting_table)
spb210_plot <- ggplot(
  spb210_plotting_table,
  aes(x = x, y = y)
) +
  geom_point() +
  geom_errorbar(
    aes(
      xmin = x - x_uncertainty,
      xmax = x + x_uncertainty
    ),
    width = 0.2
  ) +
  labs(
    title = "Supported Pb-210 vs Depth",
    x = "Supported Pb-210 (dpm/g)",
    y = "Mid Depth (cm)"
  ) +
  dynamic_scale(
    spb210_plotting_table,
    "x"
  ) +
  dynamic_scale(
    spb210_plotting_table,
    "y",
    axis = "y",
    reverse = TRUE
  ) +
  theme_minimal()

print(spb210_plot)
dev.off()
# total pb-210
pdf(file.path(folder_name, "total_pb210.pdf"), width = 8, height = 6)
tpb210_plotting_table=array(data=NA, dim=c(nrow(mass_table)/2,3), dimnames=list(NULL, c("x","y","x_uncertainty")))
tpb210_plotting_table[,1]=concentrations_table[seq(2, nrow(concentrations_table), by=2),2]
tpb210_plotting_table[,2]=concentrations_table[seq(2, nrow(concentrations_table), by=2),1]
tpb210_plotting_table[,3]=concentrations_table[seq(2, nrow(concentrations_table), by=2),3]
tpb210_plotting_table=as.data.frame(tpb210_plotting_table)
tpb210_plot <- ggplot(
  tpb210_plotting_table,
  aes(x = x, y = y)
) +
  geom_point() +
  geom_errorbar(
    aes(
      xmin = x - x_uncertainty,
      xmax = x + x_uncertainty
    ),
    width = 0.2
  ) +
  labs(
    title = "Total Pb-210 vs Depth",
    x = "Total Pb-210 (dpm/g)",
    y = "Mid Depth (cm)"
  ) +
  dynamic_scale(
    tpb210_plotting_table,
    "x"
  ) +
  dynamic_scale(
    tpb210_plotting_table,
    "y",
    axis = "y",
    reverse = TRUE
  ) +
  theme_minimal()

print(tpb210_plot)
dev.off()
# total log scale pb-210
pdf(file.path(folder_name, "excess_log_pb210.pdf"), width = 8, height = 6)
tlogpb210_plotting_table=array(data=NA, dim=c(nrow(mass_table)/2,3), dimnames=list(NULL, c("x","y","x_uncertainty")))
tlogpb210_plotting_table[,1]=concentrations_table[seq(2, nrow(concentrations_table), by=2),13]
tlogpb210_plotting_table[,2]=concentrations_table[seq(2, nrow(concentrations_table), by=2),1]
tlogpb210_plotting_table[,3]=concentrations_table[seq(2, nrow(concentrations_table), by=2),14]
tlogpb210_plotting_table=as.data.frame(tlogpb210_plotting_table)
tlogpb210_plot <- ggplot(
  tlogpb210_plotting_table,
  aes(x = x, y = y)
) +
  geom_point() +
  geom_errorbar(
    aes(
      xmin = x - x_uncertainty,
      xmax = x + x_uncertainty
    ),
    width = 0.2
  ) +
  labs(
    title = "Excess Log-scale Pb-210 vs Depth",
    x = "Excess Log-scale Pb-210 (dpm/g)",
    y = "Mid Depth (cm)"
  ) +
  dynamic_scale(
    tlogpb210_plotting_table,
    "x"
  ) +
  dynamic_scale(
    tlogpb210_plotting_table,
    "y",
    axis = "y",
    reverse = TRUE
  ) +
  theme_minimal()

print(tlogpb210_plot)
dev.off()

# mass depth version of above plots
# excess pb-210 mass depth
pdf(file.path(folder_name, "excess_pb210_md.pdf"), width = 8, height = 6)
epb210_md_plotting_table=array(data=NA, dim=c(nrow(mass_table)/2,3), dimnames=list(NULL, c("x","y","x_uncertainty")))
epb210_md_plotting_table[,1]=na.omit(concentrations_table[,6])
epb210_md_plotting_table[,2]=na.omit(concentrations_table[,12])
epb210_md_plotting_table[,3]=na.omit(concentrations_table[,7])
epb210_md_plotting_table=as.data.frame(epb210_md_plotting_table)
epb210_md_plot <- ggplot(
  epb210_md_plotting_table,
  aes(x = x, y = y)
) +
  geom_point() +
  geom_errorbar(
    aes(
      xmin = x - x_uncertainty,
      xmax = x + x_uncertainty
    ),
    width = 0.2
  ) +
  labs(
    title = "Excess Pb-210 vs Cumulative Mass Depth",
    x = "Excess Pb-210 (dpm/g)",
    y = "Cumulative Mass Depth (g/cm^2)"
  ) +
  dynamic_scale(
    epb210_md_plotting_table,
    "x"
  ) +
  dynamic_scale(
    epb210_md_plotting_table,
    "y",
    axis = "y",
    reverse = TRUE
  ) +
  theme_minimal()

print(epb210_md_plot)
dev.off()
# supported pb-210 mass depth
pdf(file.path(folder_name, "supported_pb210_md.pdf"), width = 8, height = 6)
spb210_md_plotting_table=array(data=NA, dim=c(nrow(mass_table)/2,3), dimnames=list(NULL, c("x","y","x_uncertainty")))
spb210_md_plotting_table[,1]=concentrations_table[seq(2, nrow(concentrations_table), by=2),4]
spb210_md_plotting_table[,2]=concentrations_table[seq(2, nrow(concentrations_table), by=2),12]
spb210_md_plotting_table[,3]=concentrations_table[seq(2, nrow(concentrations_table), by=2),5]
spb210_md_plotting_table=as.data.frame(spb210_md_plotting_table)
spb210_md_plot <- ggplot(
  spb210_md_plotting_table,
  aes(x = x, y = y)
) +
  geom_point() +
  geom_errorbar(
    aes(
      xmin = x - x_uncertainty,
      xmax = x + x_uncertainty
    ),
    width = 0.2
  ) +
  labs(
    title = "Supported Pb-210 vs Cumulative Mass Depth",
    x = "Supported Pb-210 (dpm/g)",
    y = "Cumulative Mass Depth (g/cm^2)"
  ) +
  dynamic_scale(
    spb210_md_plotting_table,
    "x"
  ) +
  dynamic_scale(
    spb210_md_plotting_table,
    "y",
    axis = "y",
    reverse = TRUE
  ) +
  theme_minimal()

print(spb210_md_plot)
dev.off()
# total pb-210 mass depth 
pdf(file.path(folder_name, "total_pb210_md.pdf"), width = 8, height = 6)
tpb210_md_plotting_table=array(data=NA, dim=c(nrow(mass_table)/2,3), dimnames=list(NULL, c("x","y","x_uncertainty")))
tpb210_md_plotting_table[,1]=concentrations_table[seq(2, nrow(concentrations_table), by=2),2]
tpb210_md_plotting_table[,2]=concentrations_table[seq(2, nrow(concentrations_table), by=2),12]
tpb210_md_plotting_table[,3]=concentrations_table[seq(2, nrow(concentrations_table), by=2),3]
tpb210_md_plotting_table=as.data.frame(tpb210_md_plotting_table)
tpb210_md_plot <- ggplot(
  tpb210_md_plotting_table,
  aes(x = x, y = y)
) +
  geom_point() +
  geom_errorbar(
    aes(
      xmin = x - x_uncertainty,
      xmax = x + x_uncertainty
    ),
    width = 0.2
  ) +
  labs(
    title = "Total Pb-210 vs Cumulative Mass Depth",
    x = "Total Pb-210 (dpm/g)",
    y = "Cumulative Mass Depth (g/cm^2)"
  ) +
  dynamic_scale(
    tpb210_md_plotting_table,
    "x"
  ) +
  dynamic_scale(
    tpb210_md_plotting_table,
    "y",
    axis = "y",
    reverse = TRUE
  ) +
  theme_minimal()

print(tpb210_md_plot)
dev.off()
# Excess log scale pb-210 mass depth
pdf(file.path(folder_name, "excess_log_pb210_md.pdf"), width = 8, height = 6)
tlogpb210_md_plotting_table=array(data=NA, dim=c(nrow(mass_table)/2,3), dimnames=list(NULL, c("x","y","x_uncertainty")))
tlogpb210_md_plotting_table[,1]=log(na.omit(concentrations_table[,6]))
tlogpb210_md_plotting_table[,2]=na.omit(concentrations_table[,12])
tlogpb210_md_plotting_table[,3]=na.omit(concentrations_table[,7]/concentrations_table[,6])
tlogpb210_md_plotting_table=as.data.frame(tlogpb210_md_plotting_table)
tlogpb210_md_plot <- ggplot(
  tlogpb210_md_plotting_table,
  aes(x = x, y = y)
) +
  geom_point() +
  geom_errorbar(
    aes(
      xmin = x - x_uncertainty,
      xmax = x + x_uncertainty
    ),
    width = 0.2
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    color = "red"
  ) +
  labs(
    title = "Excess Log-scale Pb-210 vs Cumulative Mass Depth",
    x = "Excess Log-scale Pb-210 (dpm/g)",
    y = "Cumulative Mass Depth (g/cm^2)"
  ) +
  dynamic_scale(
    tlogpb210_md_plotting_table,
    "x"
  ) +
  dynamic_scale(
    tlogpb210_md_plotting_table,
    "y",
    axis = "y",
    reverse = TRUE
  ) +
  theme_minimal()

print(tlogpb210_md_plot)
dev.off()
#make constants table to save
constants_table <- data.frame(
  constant = c(
    "day",
    "month",
    "year",
    "full numeric year",
    "core diameter (cm)",
    "core diameter uncertainty (cm)",
    "surface area (cm^2)",
    "surface area uncertainty (cm^2)",
    "interval thickness (cm)",
    "full core length (cm)",
    "Pb-210 half-life (yr)",
    "Pb-210 half-life uncertainty (yr)",
    "Pb-210 disintegration constant (yr^-1)",
    "Pb-210 disinetragtion constant uncertainty (yr^-1)"
  ),
  value = I(vector("list", 14)),   # list-column
  stringsAsFactors = FALSE
)

# 2) Fill in scalar entries as 1-length numerics
constants_table$value[[1]]  <- day
constants_table$value[[2]]  <- month
constants_table$value[[3]]  <- year
constants_table$value[[4]]  <- numeric_date
constants_table$value[[5]]  <- internal_diameter
constants_table$value[[6]]  <- internal_diameter_uncer
constants_table$value[[7]]  <- surface_area
constants_table$value[[8]]  <- surface_area_uncer

# 3) This row can now hold a numeric vector (float list)
if (exists("interval_thicknesses_list")) {
  constants_table$value[[9]]  <- interval_thicknesses_list   # vector allowed here
} else {
  constants_table$value[[9]]  <- interval_thickness
}

constants_table$value[[10]] <- core_length
constants_table$value[[11]] <- half_life
constants_table$value[[12]] <- half_life_uncer
constants_table$value[[13]] <- decay_const
constants_table$value[[14]] <- decay_const_uncer

# 4) Make a CSV-friendly copy (so write.csv still works cleanly)
constants_table_csv <- constants_table

constants_table_csv$value <- vapply(
  constants_table$value,
  function(x) {
    # x might be a numeric, or a list element containing a numeric/vector
    if (is.list(x)) {
      x <- unlist(x, recursive = FALSE, use.names = FALSE)
    }
    # Now x is atomic: length 1 (scalar) or >1 (vector)
    if (length(x) <= 1) {
      as.character(x)
    } else {
      paste(x, collapse = ";")  # works for either case
    }
  },
  character(1)
)

constants_file_path <- file.path(folder_name, "constants_table.xlsx")
write_xlsx(constants_table_csv, constants_file_path)
#start of CFCS section
CFCS_table=array(data=NA, dim=c(nrow(concentrations_table)/2,6), dimnames=list(NULL,c("zi","mi","ln(Ci)","Time yr","Date (CE)", "Time Uncertainty (yr)")))
colnames(CFCS_table) <- gsub("delta", "\u0394", colnames(CFCS_table))
CFCS_table[,1]=na.omit(concentrations_table[,1])
CFCS_table[,2]=na.omit(concentrations_table[,12])
for(i in 1:nrow(concentrations_table)){
  if (i%%2==0){
    CFCS_table[(i/2),3]=concentrations_table[i,13]
  }
}
CFCS_table=as.data.frame(CFCS_table)

#generate compaction corrected depths for alternative SAR calculations
zi_correction <- readline(prompt ="Do you want to calculate compaction corrected linear depths and dry bulk densities for SAR calculations? (TRUE or FALSE)? ")
zi_correction <- tolower(trimws(zi_correction))
if (zi_correction %in% c("true", "t", "yes", "y")) {
  print("select porosity data")
  file_path_porosity <- file.choose()
  porosity_data <- read_excel(file_path_porosity,col_names = TRUE,col_types = "numeric", .name_repair = "unique_quiet")
  porosity_data <- as.data.frame(porosity_data)
  measured_phi <- as.numeric(porosity_data[[1]])
  measured_phi_sd <- as.numeric(porosity_data[[2]])
  n_measured <- length(measured_phi)
  zi_correction_table <-as.data.frame(mass_table[, 1:3])
  zi_correction_table$`Porosity` <- NA_real_
  zi_correction_table$`Porosity Uncertainty` <- NA_real_
  zi_correction_table$`Uncompacted Section Thickness (cm)` <-NA_real_
  zi_correction_table$`Uncompacted Section Thickness Uncertainty (cm)` <-NA_real_
  zi_correction_table$`Uncompacted Top of Interval (cm)` <-NA_real_
  zi_correction_table$`Uncompacted Top of Interval Uncertainty (cm)` <-NA_real_
  zi_correction_table$`Uncompacted Mid-Depth (cm)` <-NA_real_
  zi_correction_table$`Uncompacted Mid-Depth Uncertainty (cm)` <-NA_real_
  
  all_rows <- seq_len(nrow(zi_correction_table))
  section_rows <- all_rows[all_rows %% 2 == 0]
  top_rows <- all_rows[all_rows %% 2 == 1]
  n_sections <- length(section_rows)
  
  # Assign measured porosity values to corresponding sections
  n_insert <- min(n_sections, n_measured)
  
  if (n_insert > 0) {
    measured_sections <- seq_len(n_insert)
    
    zi_correction_table[section_rows[measured_sections], "Porosity"] <-
      measured_phi[measured_sections]
    
    zi_correction_table[section_rows[measured_sections], "Porosity Uncertainty"] <-
      measured_phi_sd[measured_sections]
  }
  measured_porosity=readline(prompt="Do you have a measured surface layer porosity value? (TRUE or FALSE)")
  
  if (measured_porosity==TRUE){
    
    # PATHWAY 1: Independently measured surface porosity
    
    phi_0=as.numeric(readline(prompt="Enter surface layer porosity"))
    phi_0_uncer=as.numeric(readline(prompt="Enter surface layer porosity uncertainty"))
    
    phi_0_variance <- phi_0_uncer^2
    Cov_phi0_measured <- numeric(n_measured)
    
  }else{
    
    # PATHWAY 2: Regression-derived surface porosity
    
    porosity_table=array(data=NA, dim=c(10,2), dimnames=list(NULL, c("zi","phi(i)")))
    
    porosity_table <- data.frame(
      zi = as.numeric(zi_correction_table[seq(2, 20, by = 2), 2]),
      phi_i = as.numeric(zi_correction_table[seq(2, 20, by = 2), "Porosity"])
    )
    
    porosity_model_list <- list()
    model_summaries_porosity <- list()
    model_source_rows <- list()
    
    for (extent in 3:10) {
      
      candidate <- porosity_table[seq_len(extent), ]
      
      valid <- is.finite(candidate[, "zi"]) &
        is.finite(candidate[, "phi_i"])
      
      fit_data <- candidate[valid, , drop = FALSE]
      
      result <- data.frame(
        Model = paste0("Indices 1-", extent),
        N = nrow(fit_data),
        Used_indices = paste(which(valid), collapse = ", "),
        Intercept = NA_real_,
        Slope = NA_real_,
        Intercept_SE = NA_real_,
        Slope_SE = NA_real_,
        R.squared = NA_real_,
        Adj.R.squared = NA_real_,
        Status = "Insufficient valid observations",
        stringsAsFactors = FALSE
      )
      
      if (nrow(fit_data) >= 3 &&
          length(unique(fit_data$zi)) >= 2) {
        
        model <- lm(phi_i ~ zi, data = fit_data)
        
        porosity_model_list[[as.character(extent - 2)]] <- model
        
        model_source_rows[[as.character(extent - 2)]] <-
          section_rows[which(valid)]
        
        sm <- summary(model)
        coefs <- sm$coefficients
        
        result$Intercept <- unname(coef(model)[1])
        result$Slope <- unname(coef(model)[2])
        result$Intercept_SE <- coefs[1, "Std. Error"]
        result$Slope_SE <- coefs[2, "Std. Error"]
        result$R.squared <- sm$r.squared
        result$Adj.R.squared <- sm$adj.r.squared
        result$Status <- "Fitted"
      }
      
      model_summaries_porosity[[as.character(extent)]] <- result
    }
    
    model_summary_table <- do.call(rbind, model_summaries_porosity)
    model_summaries_porosity <- model_summary_table
    rownames(model_summaries_porosity) <- NULL
    
    print("Determine surface layer porosity")
    print(porosity_table)
    print(model_summaries_porosity[,-(2:3)])
    
    selection_complete <- FALSE
    
    while (selection_complete == FALSE) {
      
      best_model_index <- readline(prompt = "Select model (1-8) ")
      
      best_model <- porosity_model_list[[best_model_index]]
      
      print("Estimated surface porosity")
      print(unname(coef(best_model)[1]))
      
      selection_complete <- readline(
        prompt = "Accept model selection? (TRUE or FALSE) "
      )
    }
    
    best_model <- porosity_model_list[[best_model_index]]
    best_source_rows <- model_source_rows[[best_model_index]]
    
    beta <- coef(best_model)
    phi_0 <- unname(beta["(Intercept)"])
  }
  
  
  # COMMON PATHWAY: Down-core porosity and covariance
  
  if (length(measured_phi_sd) != n_measured) {
    stop("Porosity measurements and uncertainties have different lengths.")
  }
  
  valid_measurement <- is.finite(measured_phi)
  
  measured_phi_for_math <- measured_phi
  measured_phi_for_math[!valid_measurement] <- 0
  
  if (any(!is.finite(measured_phi_sd[valid_measurement])) ||
      any(measured_phi_sd[valid_measurement] < 0)) {
    stop("Invalid uncertainty for a measured porosity.")
  }
  
  measurement_variance <- numeric(n_measured)
  
  measurement_variance[valid_measurement] <-
    measured_phi_sd[valid_measurement]^2
  
  Sigma_measured <- diag(
    measurement_variance,
    nrow = n_measured,
    ncol = n_measured
  )
  
  # Construct measurement weight matrix
  
  W <- matrix(0, nrow = n_sections, ncol = n_measured)
  
  n_direct <- min(n_sections, n_measured)
  
  if (n_direct > 0) {
    for (i in seq_len(n_direct)) {
      if (valid_measurement[i]) {
        W[i, i] <- 1
      }
    }
  }
  
  W_original <- W
  
  # Interpolate missing porosity using original measurements only
  
  for (i in seq_len(n_sections)) {
    
    if (all(W_original[i, ] == 0)) {
      
      j <- i - 1
      
      while (j >= 1 && all(W_original[j, ] == 0)) {
        j <- j - 1
      }
      
      k <- i + 1
      
      while (k <= n_sections && all(W_original[k, ] == 0)) {
        k <- k + 1
      }
      
      if (j >= 1 && k <= n_sections) {
        
        depth_i <- zi_correction_table[
          section_rows[i], "Mid Depth zi cm"
        ]
        
        depth_j <- zi_correction_table[
          section_rows[j], "Mid Depth zi cm"
        ]
        
        depth_k <- zi_correction_table[
          section_rows[k], "Mid Depth zi cm"
        ]
        
        denominator <- depth_k - depth_j
        
        if (!all(is.finite(c(depth_i, depth_j, depth_k))) ||
            denominator == 0) {
          stop(paste("Invalid interpolation depths at section", i))
        }
        
        fraction <- (depth_i - depth_j) / denominator
        
        if (fraction < 0 || fraction > 1) {
          stop(paste("Invalid interpolation fraction at section", i))
        }
        
        W[i, ] <-
          (1 - fraction) * W_original[j, ] +
          fraction * W_original[k, ]
      }
    }
  }
  
  if (any(rowSums(abs(W)) == 0)) {
    stop("Some sections have no valid porosity measurements bounding them.")
  }
  
  # Down-core porosity vector and covariance
  
  phi <- as.numeric(W %*% measured_phi_for_math)
  
  Sigma_phi <- W %*% Sigma_measured %*% t(W)
  Sigma_phi <- (Sigma_phi + t(Sigma_phi)) / 2
  
  phi_uncertainty <- sqrt(pmax(diag(Sigma_phi), 0))
  
  zi_correction_table[section_rows, "Porosity"] <- phi
  zi_correction_table[section_rows, "Porosity Uncertainty"] <- phi_uncertainty
  
  
  # SURFACE POROSITY COVARIANCE
  
  if (measured_porosity != TRUE) {
    
    model_x <- model.matrix(best_model)
    
    regression_rows <- best_source_rows
    
    measured_section_lookup <- rep(NA_integer_, n_sections)
    
    if (n_direct > 0) {
      for (i in seq_len(n_direct)) {
        if (valid_measurement[i]) {
          measured_section_lookup[i] <- i
        }
      }
    }
    
    regression_measurement_index <- vapply(
      regression_rows,
      function(table_row) {
        
        section_index <- match(table_row, section_rows)
        
        if (is.na(section_index)) return(NA_integer_)
        
        measured_section_lookup[section_index]
      },
      integer(1)
    )
    
    if (anyNA(regression_measurement_index) ||
        any(regression_measurement_index < 1) ||
        any(regression_measurement_index > n_measured)) {
      stop("Regression rows do not map to valid measured porosity indices.")
    }
    
    Sigma_regression_y <- diag(
      measured_phi_sd[regression_measurement_index]^2
    )
    
    # OLS coefficient weights
    
    A_ols <- solve(t(model_x) %*% model_x) %*% t(model_x)
    
    # Regression coefficient covariance
    
    Sigma_beta <- A_ols %*%
      Sigma_regression_y %*%
      t(A_ols)
    
    phi_0_variance <- Sigma_beta[1, 1]
    phi_0_uncer <- sqrt(phi_0_variance)
    
    # Covariance between surface porosity and original measurements
    
    intercept_weights <- A_ols[1, ]
    
    Cov_phi0_regression_y <- as.numeric(
      intercept_weights %*% Sigma_regression_y
    )
    
    Cov_phi0_measured <- numeric(n_measured)
    
    Cov_phi0_measured[regression_measurement_index] <-
      Cov_phi0_regression_y
  }
  
  # Propagate surface porosity covariance through W
  
  Cov_phi0_phi <- as.numeric(W %*% Cov_phi0_measured)
  
  
  # FULL INPUT COVARIANCE MATRIX
  
  Sigma_input <- matrix(
    0,
    nrow = n_sections + 1,
    ncol = n_sections + 1
  )
  
  Sigma_input[
    seq_len(n_sections),
    seq_len(n_sections)
  ] <- Sigma_phi
  
  Sigma_input[
    seq_len(n_sections),
    n_sections + 1
  ] <- Cov_phi0_phi
  
  Sigma_input[
    n_sections + 1,
    seq_len(n_sections)
  ] <- Cov_phi0_phi
  
  Sigma_input[
    n_sections + 1,
    n_sections + 1
  ] <- phi_0_variance
  
  Sigma_input <- (Sigma_input + t(Sigma_input)) / 2
  #decompaction
  # h'_i =
  #
  #       h_i (1 - phi_i)
  #       -----------------
  #          (1 - phi_0)
  h <- as.numeric(zi_correction_table[section_rows,3])
  # Corrected thickness
  h_prime <-h *(1 - phi) /(1 - phi_0)
  # Store
  zi_correction_table[section_rows,"Uncompacted Section Thickness (cm)"] <-h_prime
  # JACOBIAN OF CORRECTED THICKNESS
  # Since h is exact:
  #
  # d h'_i / d phi_i
  #
  #       = -h_i / (1 - phi_0)
  # d h'_i / d phi_0
  #
  #       = h'_i / (1 - phi_0)
  J_hprime <- matrix(0,nrow = n_sections,ncol = n_sections + 1)
  for (i in seq_len(n_sections)) {
    J_hprime[i,i] <- -h[i] /(1 - phi_0)
  }
  # Every section depends on the common phi_0
  J_hprime[,n_sections + 1] <- h_prime /(1 - phi_0)
  #covariance of decompacted section thickness
  Sigma_hprime <-J_hprime %*% Sigma_input %*% t(J_hprime)
  Sigma_hprime <- (Sigma_hprime + t(Sigma_hprime)) / 2
  #section thickness uncertainty
  h_prime_uncertainty <- sqrt(pmax(diag(Sigma_hprime),0))
  zi_correction_table[section_rows,"Uncompacted Section Thickness Uncertainty (cm)"] <-h_prime_uncertainty
  #decompacted layer depths
  z_prime <- numeric(n_sections)
  z_prime[1] <- 0
  if (n_sections > 1) {
    z_prime[2:n_sections] <-cumsum(h_prime)[1:(n_sections - 1)]
  }
  # Assign top-of-interval depths to the corresponding boundary rows
  zi_correction_table[top_rows[seq_len(n_sections)],
                      "Uncompacted Top of Interval (cm)"] <- z_prime
  # layer depths covariance
  L <- matrix(0,nrow = n_sections,ncol = n_sections)
  if (n_sections > 1) {
    for (i in 2:n_sections) {
      L[i,1:(i - 1)] <- 1
    }
  }
  
  Sigma_top <-L %*% Sigma_hprime %*% t(L)
  Sigma_top <- (Sigma_top + t(Sigma_top)) / 2
  #uncompacted layer depth uncertainty
  z_prime_uncertainty <-sqrt(pmax(diag(Sigma_top),0))
  zi_correction_table[top_rows[seq_len(n_sections)],
                      "Uncompacted Top of Interval Uncertainty (cm)"] <- z_prime_uncertainty
  #uncompacted midpoint depths
  mid_prime <-z_prime +0.5 * h_prime
  zi_correction_table[section_rows,"Uncompacted Mid-Depth (cm)"] <-mid_prime
  #covariance of midpoints
  # m' = z' + 0.5 h'
  #       = M h'
  M <- L
  for (i in seq_len(n_sections)) {
    M[i,i] <- 0.5
  }
  Sigma_mid <- M %*% Sigma_hprime %*% t(M)
  Sigma_mid <- (Sigma_mid + t(Sigma_mid)) / 2
  #midpoint depth uncertainty
  mid_prime_uncertainty <- sqrt(pmax(diag(Sigma_mid),0))
  zi_correction_table[section_rows,"Uncompacted Mid-Depth Uncertainty (cm)"] <- mid_prime_uncertainty
  # SAVE COVARIANCE MATRICES AS ATTRIBUTES
  # These are useful if the corrected depths are subsequently
  # used in another uncertainty calculation.
  attr(zi_correction_table,"Surface_Porosity") <- phi_0
  attr(zi_correction_table,"Surface_Porosity_Uncertainty") <- phi_0_uncer
  attr(zi_correction_table,"Cov_phi0_phi") <- Cov_phi0_phi
  attr(zi_correction_table,"Porosity_Covariance") <- Sigma_phi
  attr(zi_correction_table,"Corrected_Thickness_Covariance") <- Sigma_hprime
  attr(zi_correction_table,"Corrected_TopDepth_Covariance") <- Sigma_top
  attr(zi_correction_table,"Corrected_MidDepth_Covariance") <- Sigma_mid
  folder_name_depth = "Depth Correction"
  dir.create(folder_name_depth)
  depth_corr_path <- file.path(folder_name_depth, "uncompacted_depth_table.xlsx")
  write_xlsx(zi_correction_table, depth_corr_path)
  if (exists("model_summaries_porosity")){
  write_xlsx(model_summaries_porosity, path = "Depth Correction/porosity_model_summaries.xlsx")
  }
  }

#start of SAZ auto detection
SAZ_auto_detect=na.omit(concentrations_table[,6])
SAZ_auto_detect=as.data.frame(SAZ_auto_detect)
SAZ_auto_detect$X.2=na.omit(concentrations_table[,1])
SAZ_auto_detect=as.data.frame(SAZ_auto_detect)
# Perform analysis
detect_pb210_reworking_top <- function(x, y,
                                       hmin = 3,
                                       slope_tol = NULL,      # dpm/g per cm; if NULL, data-driven
                                       alpha_diff = 0.05,     # slope differs below vs above (want small p)
                                       min_rise = NULL) {     # require mean(below)-mean(top) >= min_rise
  
  x <- as.numeric(x); y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  
  o <- order(x)
  x <- x[o]; y <- y[o]
  
  n <- length(y)
  if (n < 2 * hmin + 1) {
    return(list(present = FALSE, idx = integer(0), reason = "too_few_points"))
  }
  
  if (is.null(slope_tol)) {
    b_global <- coef(lm(y ~ x))[2]
    noise_floor <- 0.25 * stats::mad(diff(y), constant = 1)
    slope_tol <- max(0.1 * abs(b_global), noise_floor)
    if (!is.finite(slope_tol) || slope_tol <= 0) slope_tol <- noise_floor
  }
  
  if (is.null(min_rise)) {
    min_rise <- stats::mad(y, constant = 1)
    if (!is.finite(min_rise) || min_rise < 0) min_rise <- 0
  }
  
  # j = last index of the "flat top" segment (top = 1:j, below = (j+1):n)
  j_grid <- hmin:(n - hmin)
  
  rows <- lapply(j_grid, function(j) {
    top <- 1:j
    bot <- (j + 1):n
    
    fit_top <- lm(y[top] ~ x[top])
    fit_bot <- lm(y[bot] ~ x[bot])
    
    b_top <- unname(coef(fit_top)[2])
    b_bot <- unname(coef(fit_bot)[2])
    
    se_top <- summary(fit_top)$coefficients[2, 2]
    se_bot <- summary(fit_bot)$coefficients[2, 2]
    
    df_top <- length(top) - 2
    df_bot <- length(bot) - 2
    
    # slope difference test (want significant)
    tdiff <- (b_bot - b_top) / sqrt(se_top^2 + se_bot^2)
    df_diff <- (se_top^2 + se_bot^2)^2 / ((se_top^4 / df_top) + (se_bot^4 / df_bot))
    p_diff <- 2 * pt(abs(tdiff), df = df_diff, lower.tail = FALSE)
    
    mu_top <- mean(y[top])
    mu_bot <- mean(y[bot])
    rise <- mu_bot - mu_top
    
    data.frame(
      j = j,
      x_end_top = x[j],
      n_top = length(top),
      n_bot = length(bot),
      slope_top = b_top,
      abs_slope_top = abs(b_top),
      slope_bot = b_bot,
      p_diff = p_diff,
      rise = rise
    )
  })
  
  tab <- do.call(rbind, rows)
  
  keep <- with(tab,
               abs_slope_top <= slope_tol &
                 p_diff <= alpha_diff &
                 rise >= min_rise)
  
  if (!any(keep)) {
    return(list(
      present = FALSE,
      idx = integer(0),
      slope_tol = slope_tol,
      min_rise = min_rise,
      candidates = tab
    ))
  }
  
  cand <- tab[keep, , drop = FALSE]
  best <- cand[order(cand$j), ][1, , drop = FALSE]   # earliest end-of-flat-top
  
  list(
    present = TRUE,
    idx = best$j,              # last index in the flat top (after sorting by depth)
    x_end_top = best$x_end_top,
    slope_top = best$slope_top,
    slope_tol = slope_tol,
    min_rise = min_rise,
    best_row = best,
    candidates = tab
  )
}

# Use (y = dpm/g, x = cm)
y <- SAZ_auto_detect[, 1]
x <- SAZ_auto_detect[, 2]
res_top <- detect_pb210_reworking_top(x, y, hmin = 3)

res_top$present
res_top$idx        # end index of the flat mixed layer
res_top$x_end_top
if (res_top$present==TRUE) {
  
  suggested_SAZ_index <- res_top$idx
  cat("Suggested SAZ terminates at index:", suggested_SAZ_index, "\n")
} else {
  cat("No SAZ detected.\n")
}
# end of SAZ auto detection
CFCS_modeling_complete=FALSE
SAZ_depth=0
  surface_table=array(data=NA, dim=c(nrow(mass_table)/2,3), dimnames=list(NULL, c("x","y", "x_uncertainty")))
  surface_table[,1]=na.omit(concentrations_table[,6])
  surface_table[,2]=na.omit(concentrations_table[,1])
  surface_table[,3]=na.omit(concentrations_table[,7])
  surface_table=as.data.frame(surface_table)
  surface_table_plot <- ggplot(
    surface_table,
    aes(x = x, y = y)
  ) +
    geom_errorbar(
      aes(
        xmin = x - x_uncertainty,
        xmax = x + x_uncertainty
      ),
      width = 0.2
    ) +
    geom_point() +
    labs(
      title = "Excess Pb-210 vs Mid Depth",
      x = "Excess Pb-210 (dpm/g)",
      y = "Mid Depth (cm)"
    ) +
    dynamic_scale(
      surface_table,
      "x"
    ) +
    dynamic_scale(
      surface_table,
      "y",
      axis = "y",
      reverse = TRUE
    ) +
    theme_minimal()
  
  print(surface_table_plot)
  attributes(surface_table)$na.action <- NULL
  print("check for surface active zone")
  SAZ_determination_complete = FALSE
  
  while (!SAZ_determination_complete) {
    SAZ_view = as.numeric(readline(prompt = "How many datapoints do you wish to view for surface active zone determination? "))
    
    # Ensure SAZ_view is valid
    if (SAZ_view <= 0 || SAZ_view > nrow(surface_table)) {
      cat("Please enter a valid number of datapoints.\n")
      next
    }
    
    # Create a temporary table for viewing
    temp_surface_table = surface_table[1:SAZ_view, 1:2]
    rownames(temp_surface_table) = as.character(seq(1, SAZ_view))
    colnames(temp_surface_table) = c("Excess Pb-210 (dpm/g)", "Mid Depth (cm)")
    attributes(temp_surface_table)$na.action <- NULL
    print(temp_surface_table)
    
    surface_active_zone_check = readline(prompt = "Does this core have a surface active zone? (True or False) ")
    
    if (tolower(surface_active_zone_check) == "false") {
      SAZ_determination_complete = TRUE
    } else if (tolower(surface_active_zone_check) == "true") {
      SAZ_depth = as.numeric(readline(prompt = "At what index does the surface active zone terminate? "))
      
      # Ensure SAZ_depth is valid
      if (SAZ_depth < 1 || SAZ_depth > nrow(CFCS_table)) {
        cat("Please enter a valid index for SAZ depth.\n")
        next
      }
      
      accept_determination = readline(prompt = "Accept current surface active zone determination? (True or False) ")
      
      if (tolower(accept_determination) == "true") {
        # Only trim the CFCS_table after acceptance
        CFCS_table = CFCS_table[SAZ_depth:nrow(CFCS_table), ]
        SAZ_determination_complete = TRUE  # Exit the loop after accepting the determination
      }
    }
  }
print("start of CFCS analysis")
print("MAR modeling")
while (CFCS_modeling_complete==FALSE){
  MAR_model_plot <- ggplot(
    CFCS_table,
    aes(
      x = .data[[names(CFCS_table)[2]]],
      y = .data[[names(CFCS_table)[3]]]
    )
  ) +
    geom_point(
      color = "blue",
      size = 2
    ) +
    labs(
      x = "Cumulative Mass Depth (g/cm^2)",
      y = "log-scale Excess Pb-210 Activity (dpm/g)",
      title = "MAR Modeling"
    ) +
    dynamic_scale(
      CFCS_table,
      names(CFCS_table)[2]
    ) +
    dynamic_scale(
      CFCS_table,
      names(CFCS_table)[3],
      axis = "y"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  
  print(MAR_model_plot)
  changepoint_detection_model = lm(`ln(Ci)` ~ mi, data = CFCS_table, na.action = na.exclude)
  number_of_models=as.numeric(readline(prompt = "How many models do you wish to try? "))
  while(number_of_models>=nrow(CFCS_table)/3){
    number_of_models=as.numeric(readline(prompt = "Too many models suggested! Enter lower number "))
  }
  number_of_changepoints=number_of_models-1
  changepoint_list=c()
  if (number_of_changepoints>0){
    for (i in 1:number_of_changepoints){
      changepoint_list[i]=as.numeric(readline(prompt = "Enter cumulative mass depth change-point(s) "))
    }
    if (any(changepoint_list > CFCS_table[nrow(CFCS_table), 2])) {
      print("Entered value(s) greater than maximum cumulative mass depth!")
      changepoint_list[i]=as.numeric(readline(prompt = "Enter cumulative mass depth change-point(s) "))
    }
    segmented_model <- segmented(changepoint_detection_model, seg.Z = ~ mi, psi = list(mi = changepoint_list))
    estimated_breakpoints <- summary(segmented_model)$psi[, "Est."]
    estimated_breakpoints=unlist(estimated_breakpoints, use.names = FALSE)
    MAR_indexes <- numeric(length(estimated_breakpoints))
    for (i in 1:length(estimated_breakpoints)) {
      
      MAR_indexes[i]<- which.min(abs(CFCS_table[,2]-estimated_breakpoints[i]))
    }
    print("estimated change point(s)")
    print(estimated_breakpoints)
    extents=c(1)
    extents <- c(extents, lapply(MAR_indexes, function(x) x))
    extents[length(extents)+1]=nrow(CFCS_table)
    print("nearest model data values")
    model_data_vals=c()
    for (i in 1:length(extents)){
      model_data_vals[[i]]=CFCS_table[extents[[i]],2]
    }
    cat(unlist(model_data_vals), sep = " ", fill=TRUE)
    print("model data indexes")
    list_of_models=c()
    extents <- unlist(extents, use.names = FALSE)
    print(extents)
    for (i in 1:number_of_models) {
      
      list_of_models[[i]] <- lm(
        CFCS_table[extents[i]:extents[i + 1], 3] ~
          CFCS_table[extents[i]:extents[i + 1], 2],
        na.action = na.exclude
      )
      
      print(paste("Model Adj. R square", i))
      print(summary(list_of_models[[i]])$adj.r.squared)
      
      print(paste("Model MAR (g/(cm^2 y))", i))
      print(-decay_const / summary(list_of_models[[i]])$coefficients[2, 1])
      
      if (
        (-decay_const /
         summary(list_of_models[[i]])$coefficients[2, 1]) < 0
      ) {
        print(
          "model assumptions violated! Negative MAR generated, choose different MAR model fit"
        )
      }
    }
    
    # Create data frame containing fitted lines
    lines_df <- data.frame()
    
    for (i in seq_along(list_of_models)) {
      
      model <- list_of_models[[i]]
      
      # Retrieve the full x-range for this model segment
      x_vals <- CFCS_table[extents[i]:extents[i + 1], 2]
      
      # Extract the original regression coefficients
      coefficients <- coef(model)
      
      # Calculate predicted y-values across the entire x-range
      y_vals <- coefficients[1] + coefficients[2] * x_vals
      
      segment_df <- data.frame(
        x = x_vals,
        y = y_vals,
        group = as.factor(i)
      )
      
      # Retain all valid x-coordinates
      segment_df <- segment_df[is.finite(segment_df$x), ]
      
      lines_df <- rbind(lines_df, segment_df)
    }
    
    # Create final MAR plot
    MAR_model_plot <- ggplot(
      CFCS_table,
      aes(x = mi, y = `ln(Ci)`)
    ) +
      geom_point(
        color = "blue"
      ) +
      geom_line(
        data = lines_df,
        aes(
          x = x,
          y = y,
          group = group,
          color = group
        ),
        linewidth = 1.2,
        show.legend = TRUE
      ) +
      labs(
        title = "MAR Modeling",
        x = "Cumulative Mass Depth (g/cm^2)",
        y = "log-scale Excess Pb-210 Activity (dpm/g)",
        color = "Model"
      ) +
      dynamic_scale(
        CFCS_table,
        "mi"
      ) +
      dynamic_scale(
        CFCS_table,
        "ln(Ci)",
        axis = "y"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5)
      )
    
    print(MAR_model_plot)
  }else{
    extents=c(1,nrow(CFCS_table))
    print("single model fit")
    print(summary(changepoint_detection_model)$adj.r.squared)
    # Generate continuous fitted line
    single_line_df <- data.frame(
      mi = CFCS_table$mi,
      fitted_y = coef(changepoint_detection_model)[1] +
        coef(changepoint_detection_model)[2] * CFCS_table$mi
    )
    
    single_line_df <- single_line_df[
      is.finite(single_line_df$mi) &
        is.finite(single_line_df$fitted_y), ]
    
    single_line_df <- single_line_df[
      order(single_line_df$mi), ]
    
    print(
      MAR_model_plot +
        geom_line(
          data = single_line_df,
          aes(x = mi, y = fitted_y),
          inherit.aes = FALSE,
          color = "red",
          linewidth = 1
        )
    )
    print("single Model MAR (g/(cm^2 y))")
    print(-decay_const/summary(changepoint_detection_model)$coefficients[2, 1])
    MAR_final_slope=summary(changepoint_detection_model)$coefficients[2, "Estimate"]
    MAR_final_slope_err=summary(changepoint_detection_model)$coefficients[2, 2]
    MAR_final=-decay_const/summary(changepoint_detection_model)$coefficients[2, 1]
    MAR_final_err=(-decay_const/MAR_final_slope)*sqrt((MAR_final_slope_err/MAR_final_slope)^2+(decay_const_uncer/decay_const)^2)
    if ((-decay_const/summary(changepoint_detection_model)$coefficients[2, 1])<0){
      print("model assumptions violated! Negative MAR generated, choose different MAR model fit")
    }
  }
  CFCS_manual_eval=readline(prompt = "Do you wish to manually set CFCS models? (True or False) ")
  if (CFCS_manual_eval==TRUE){
    MAR_model_plot <- ggplot(
      CFCS_table,
      aes(
        x = CFCS_table[, 2],
        y = CFCS_table[, 3]
      )
    ) +
      geom_point(
        color = "blue",
        size = 2
      ) +
      labs(
        x = "Cumulative Mass Depth (g/cm^2)",
        y = "log-scale Excess Pb-210 Activity (dpm/g)",
        title = "MAR Modeling"
      ) +
      dynamic_scale(
        CFCS_table,
        names(CFCS_table)[2]
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5)
      )
    
    print(MAR_model_plot)
    
    
    changepoint_detection_model = lm(`ln(Ci)` ~ mi, data = CFCS_table, na.action=na.exclude)
    number_of_models=as.numeric(readline(prompt = "How many models do you wish to try? "))
    while(number_of_models>=nrow(CFCS_table)/3){
      number_of_models=as.numeric(readline(prompt = "Too many models suggested! Enter lower number "))
    }
    number_of_changepoints=number_of_models-1
    changepoint_list=c()
    if (number_of_changepoints>0){
      for (i in 1:number_of_changepoints){
        changepoint_list[i]=as.numeric(readline(prompt = "Enter change-point indices "))
      }
      
    }
    extents=c(1, changepoint_list, length(CFCS_table[,1]))
    print("model data indexes")
    print(extents)
    print("model data values (g/cm^2)")
    model_data_vals <- numeric(length(extents))
    
    for (i in seq_along(extents)) {
      model_data_vals[i] <- CFCS_table[extents[i], 2]
    }
    
    print(model_data_vals)
    list_of_models=c()
    for (i in 1:number_of_models){
      list_of_models[[i]]=lm(CFCS_table[extents[i]:extents[i+1],3]~CFCS_table[extents[i]:extents[i+1],2], na.action=na.exclude)
      print(paste("Model Adj. R square",i))
      print(summary(list_of_models[[i]])$adj.r.square)
      print(paste("Model MAR (g/(cm^2 y))",i))
      print(-decay_const/summary(list_of_models[[i]])$coefficients[2, 1])
      MAR_model_plot + geom_line(aes(y = fitted(list_of_models[[i]])), color = i, linewidth = 1)
    }
    # Create data frame containing fitted lines
    lines_df <- data.frame()
    
    for (i in seq_along(list_of_models)) {
      
      model <- list_of_models[[i]]
      
      # Retrieve the full x-range for this model segment
      x_vals <- CFCS_table[extents[i]:extents[i + 1], 2]
      
      # Extract the original regression coefficients
      coefficients <- coef(model)
      
      # Calculate predicted y-values across the entire x-range
      y_vals <- coefficients[1] + coefficients[2] * x_vals
      
      segment_df <- data.frame(
        x = x_vals,
        y = y_vals,
        group = as.factor(i)
      )
      
      # Retain all valid x-coordinates
      segment_df <- segment_df[is.finite(segment_df$x), ]
      
      lines_df <- rbind(lines_df, segment_df)
    }
    
    MAR_model_plot <- ggplot(
      CFCS_table,
      aes(x = mi, y = `ln(Ci)`)
    ) +
      geom_point(
        color = "blue"
      ) +
      geom_line(
        data = lines_df,
        aes(
          x = x,
          y = y,
          group = group,
          color = group
        ),
        linewidth = 1.2,
        show.legend = TRUE
      ) +
      labs(
        title = "MAR Modeling",
        x = "Cumulative Mass Depth (g/cm^2)",
        y = "log-scale Excess Pb-210 Activity (dpm/g)",
        color = "Model"
      ) +
      dynamic_scale(
        CFCS_table,
        "mi"
      ) +
      dynamic_scale(
        CFCS_table,
        "ln(Ci)",
        axis = "y"
      ) +
      theme_minimal()
    
    print(MAR_model_plot)
  }
  CFCS_modeling_complete=readline(prompt = "Do you wish to accept current CFCS model(s)? (True or False) ")
}
print("SAR modeling")

if (surface_active_zone_check==TRUE){
  for (i in 1:length(extents)){
    extents[i]=extents[i]+SAZ_depth-1
    fixer_var=0
  }
}else{
  fixer_var=1
}
SAR_model_plot <- ggplot(
  CFCS_table,
  aes(
    x = CFCS_table[, 1],
    y = CFCS_table[, 3]
  )
) +
  geom_point(
    color = "blue",
    size = 2
  ) +
  labs(
    x = "Mid Depth (cm)",
    y = "log-scale Excess Pb-210 Activity (dpm/g)",
    title = "SAR Modeling"
  ) +
  dynamic_scale(
    CFCS_table,
    names(CFCS_table)[1]
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

print(SAR_model_plot)
changepoint_detection_model_SAR = lm(`ln(Ci)` ~ zi, data = CFCS_table, na.action=na.exclude)

if (number_of_changepoints>0){
  model_data_vals_SAR=c()
  for (i in 1:length(extents)){
    model_data_vals_SAR[[i]]=CFCS_table[(extents[[i]]-(SAZ_depth-1+fixer_var)),1]
  }
  cat(unlist(model_data_vals_SAR), sep = " ", fill=TRUE)
  print("model data indexes")
  print(extents)
  list_of_models_SAR= vector("list", number_of_models)
  for (i in 1:number_of_models){
    list_of_models_SAR[[i]]=lm(CFCS_table[(extents[i]-(SAZ_depth-1+fixer_var)):(extents[i+1]-(SAZ_depth-1+fixer_var)),3]~CFCS_table[(extents[i]-(SAZ_depth-1+fixer_var)):(extents[i+1]-(SAZ_depth-1+fixer_var)),1], na.action=na.exclude)
    print(paste("Model Adj. R square",i))
    print(summary(list_of_models_SAR[[i]])$adj.r.square)
    print(paste("Model SAR (g/(cm^2 yr)",i))
    print(-decay_const/summary(list_of_models_SAR[[i]])$coefficients[2, 1])
    SAR_model_plot + geom_line(aes(y = fitted(list_of_models_SAR[[i]])), color = i, linewidth = 1)
  }
  lines_df_SAR <- data.frame()
  
  for (i in seq_along(list_of_models_SAR)) {
    
    model_SAR <- list_of_models_SAR[[i]]
    
    x_vals_SAR <- CFCS_table[
      (extents[i]-(SAZ_depth-1+fixer_var)):
        (extents[i+1]-(SAZ_depth-1+fixer_var)),
      1
    ]
    
    coefficients <- coef(model_SAR)
    
    y_vals_SAR <- coefficients[1] +
      coefficients[2] * x_vals_SAR
    
    segment_df_SAR <- data.frame(
      x = x_vals_SAR,
      y = y_vals_SAR,
      group = as.factor(i)
    )
    
    segment_df_SAR <- segment_df_SAR[
      is.finite(segment_df_SAR$x), ]
    
    lines_df_SAR <- rbind(lines_df_SAR, segment_df_SAR)
  }
  
  SAR_model_plot <- ggplot(
    CFCS_table,
    aes(x = zi, y = `ln(Ci)`)
  ) +
    geom_point(
      color = "blue"
    ) +
    geom_line(
      data = lines_df_SAR,
      aes(
        x = x,
        y = y,
        group = group,
        color = group
      ),
      linewidth = 1.2,
      show.legend = TRUE
    ) +
    labs(
      title = "SAR Modeling",
      x = "Mid Depth (cm)",
      y = "log-scale Pb-210 Activity (dpm/g)",
      color = "Model"
    ) +
    dynamic_scale(
      CFCS_table,
      "zi"
    ) +
    theme_minimal()
  
  print(SAR_model_plot)
  
}else{
  print("single model fit")
  print(summary(changepoint_detection_model_SAR)$adj.r.squared)
  print("single model SAR (cm/yr)")
  print(-decay_const/summary(changepoint_detection_model_SAR)$coefficients[2, 1])
  # Generate continuous fitted line
  single_line_df_SAR <- data.frame(
    zi = CFCS_table$zi,
    fitted_y = coef(changepoint_detection_model_SAR)[1] +
      coef(changepoint_detection_model_SAR)[2] * CFCS_table$zi
  )
  
  single_line_df_SAR <- single_line_df_SAR[
    is.finite(single_line_df_SAR$zi) &
      is.finite(single_line_df_SAR$fitted_y), ]
  
  single_line_df_SAR <- single_line_df_SAR[
    order(single_line_df_SAR$zi), ]
  
  print(
    SAR_model_plot +
      geom_line(
        data = single_line_df_SAR,
        aes(x = zi, y = fitted_y),
        inherit.aes = FALSE,
        color = "red",
        linewidth = 1
      )
  )
  
}
if (number_of_models==1){
  extents=c(1,nrow(CFCS_table))
}

if (exists("zi_correction_table")){
  print("Uncompacted SAR modeling")
  CFCS_table$`Uncompacted Mid-Depth (cm)`=na.omit(zi_correction_table[(SAZ_depth*2):nrow(zi_correction_table),"Uncompacted Mid-Depth (cm)"])
  ccSAR_model_plot <- ggplot(
    CFCS_table,
    aes(
      x = `Uncompacted Mid-Depth (cm)`,
      y = `ln(Ci)`
    )
  ) +
    geom_point(
      color = "blue",
      size = 2
    ) +
    labs(
      x = "Uncompacted Mid Depth (cm)",
      y = "log-scale Excess Pb-210 Activity (dpm/g)",
      title = "Uncompacted ccSAR Modeling"
    ) +
    dynamic_scale(
      CFCS_table,
      "Uncompacted Mid-Depth (cm)"
    ) +
    dynamic_scale(
      CFCS_table,
      "ln(Ci)",
      axis = "y"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  
  print(ccSAR_model_plot)
  changepoint_detection_model_ccSAR = lm(`ln(Ci)` ~ `Uncompacted Mid-Depth (cm)`, data = CFCS_table, na.action=na.exclude)
  
  if (number_of_changepoints>0){
    list_of_models_ccSAR=c()
    print("model data values")
    model_data_vals_ccSAR=c()
    for (i in 1:length(extents)){
      model_data_vals_ccSAR[[i]]=CFCS_table[(extents[[i]]-(SAZ_depth-1+fixer_var)),7]
    }
    cat(unlist(model_data_vals_ccSAR), sep = " ", fill=TRUE)
    print("model data indexes")
    print(extents)
    for (i in 1:number_of_models){
      list_of_models_ccSAR[[i]]=lm(CFCS_table[(extents[i]-(SAZ_depth-1+fixer_var)):(extents[i+1]-(SAZ_depth-1+fixer_var)),3]~CFCS_table[(extents[i]-(SAZ_depth-1+fixer_var)):(extents[i+1]-(SAZ_depth-1+fixer_var)),7], na.action=na.exclude)
      print(paste("Model Adj. R square",i))
      print(summary(list_of_models_ccSAR[[i]])$adj.r.square)
      print(paste("Model ccSAR (g/(cm^2 yr)",i))
      print(-decay_const/summary(list_of_models_ccSAR[[i]])$coefficients[2, 1])
      ccSAR_model_plot + geom_line(aes(y = fitted(list_of_models_ccSAR[[i]])), color = i, linewidth = 1)
    }
    lines_df_ccSAR <- data.frame()
    
    for (i in seq_along(list_of_models_ccSAR)) {
      
      model_ccSAR <- list_of_models_ccSAR[[i]]
      
      x_vals_ccSAR <- CFCS_table[
        (extents[i]-(SAZ_depth-1+fixer_var)):
          (extents[i+1]-(SAZ_depth-1+fixer_var)),
        7
      ]
      
      coefficients <- coef(model_ccSAR)
      
      y_vals_ccSAR <- coefficients[1] +
        coefficients[2] * x_vals_ccSAR
      
      segment_df_ccSAR <- data.frame(
        x = x_vals_ccSAR,
        y = y_vals_ccSAR,
        group = as.factor(i)
      )
      
      segment_df_ccSAR <- segment_df_ccSAR[
        is.finite(segment_df_ccSAR$x) &
          is.finite(segment_df_ccSAR$y), ]
      
      lines_df_ccSAR <- rbind(
        lines_df_ccSAR,
        segment_df_ccSAR
      )
    }
    
    ccSAR_model_plot <- ggplot(
      CFCS_table,
      aes(x = `Uncompacted Mid-Depth (cm)`, y = `ln(Ci)`)
    ) +
      geom_point(
        color = "blue"
      ) +
      geom_line(
        data = lines_df_ccSAR,
        aes(
          x = x,
          y = y,
          group = group,
          color = group
        ),
        linewidth = 1.2,
        show.legend = TRUE
      ) +
      labs(
        title = " Uncompacted SAR Modeling",
        x = "Mid Depth (cm)",
        y = "log-scale Pb-210 Activity (dpm/g)",
        color = "Model"
      ) +
      dynamic_scale(
        CFCS_table,
        "Uncompacted Mid-Depth (cm)"
      ) +
      theme_minimal()
    
    print(ccSAR_model_plot)

  }else{
    print("single model fit")
    print(summary(changepoint_detection_model_ccSAR)$adj.r.squared)
    print("single model Uncompacted SAR (cm/yr)")
    print(-decay_const/summary(changepoint_detection_model_ccSAR)$coefficients[2, 1])
    # Generate continuous fitted line
    single_line_df_ccSAR <- data.frame(
      x = CFCS_table$`Uncompacted Mid-Depth (cm)`,
      fitted_y = coef(changepoint_detection_model_ccSAR)[1] +
        coef(changepoint_detection_model_ccSAR)[2] *
        CFCS_table$`Uncompacted Mid-Depth (cm)`
    )
    
    single_line_df_ccSAR <- single_line_df_ccSAR[
      is.finite(single_line_df_ccSAR$x) &
        is.finite(single_line_df_ccSAR$fitted_y), ]
    
    single_line_df_ccSAR <- single_line_df_ccSAR[
      order(single_line_df_ccSAR$x), ]
    
    print(
      ccSAR_model_plot +
        geom_line(
          data = single_line_df_ccSAR,
          aes(x = x, y = fitted_y),
          inherit.aes = FALSE,
          color = "red",
          linewidth = 1
        )
    )
    
  }
  if (number_of_models==1){
    extents=c(1,nrow(CFCS_table))
  }
}
depth_correct=0
if (surface_active_zone_check==TRUE) {
  depth_correct=SAZ_depth-1
}

model_transition_age_catch=numeric(nrow(CFCS_table))
#populate CFCS table
if (number_of_models>1){
  time_correction=0
  for (i in 1:length(list_of_models)){
    for (j in extents[i]:(extents[i+1])){
      if (time_correction==0){
        CFCS_table[(j-extents[1]+1),4]=(CFCS_table[(j-extents[1]+1),2]-CFCS_table[(extents[i]-depth_correct),2])/(-decay_const/(summary(list_of_models[[i]])$coefficients[2, "Estimate"]))+time_correction
        CFCS_table[(j-extents[1]+1),5]=numeric_date-CFCS_table[(j-extents[1]+1),4]
        model_transition_age_catch[(j-extents[1]+1)]=(CFCS_table[(j-extents[1]+1),2]-CFCS_table[(extents[i]-depth_correct),2])/(-decay_const/(summary(list_of_models[[i]])$coefficients[2, "Estimate"]))
      }else{
        CFCS_table[(j-extents[1]+1),4]=(CFCS_table[(j-extents[1]+1),2]-CFCS_table[(extents[i]-depth_correct-1),2])/(-decay_const/(summary(list_of_models[[i]])$coefficients[2, "Estimate"]))+time_correction
        CFCS_table[(j-extents[1]+1),5]=numeric_date-CFCS_table[(j-extents[1]+1),4]
        model_transition_age_catch[(j-extents[1]+1)]=(CFCS_table[(j-extents[1]+1),2]-CFCS_table[(extents[i]-depth_correct-1),2])/(-decay_const/(summary(list_of_models[[i]])$coefficients[2, "Estimate"]))
      }
      if (j>=extents[i+1]-1 && j+1!=extents[i+1]){
        time_correction=CFCS_table[(j-extents[1]),4]
      }
    }
  }
}else{
  for (i in extents[1]:extents[2]){
    CFCS_table[(i-extents[1]+1),4]=(CFCS_table[(i-extents[1]+1),2]-CFCS_table[1,2])/(-decay_const/(summary(changepoint_detection_model)$coefficients[2, "Estimate"]))
    CFCS_table[(i-extents[1]+1),5]=numeric_date-CFCS_table[(i-extents[1]+1),4]
  }
}

# calculate the age uncertainty
#find the age uncertainty for each overlying model
all_transition_age_unc_cont=numeric(number_of_changepoints)
if(exists("list_of_models")==FALSE){
  to_model_num=1
}else{
  to_model_num=length(list_of_models)
}
for (i in 1:to_model_num){
  for (j in extents[i]:extents[i+1] )
    if (j %in% (extents[-c(1, length(extents))] - 1)) {
      all_transition_age_unc_cont[i]=model_transition_age_catch[(j-extents[1]+1)]
      all_transition_age_unc_cont[i]=all_transition_age_unc_cont[i]*((-decay_const/summary(list_of_models[[i]])$coefficients[2, 1])*sqrt(((-summary(list_of_models[[i]])$coefficients[2, 2])/(summary(list_of_models[[i]])$coefficients[2,1]))^2+(decay_const_uncer/decay_const)^2)/(-decay_const/summary(list_of_models[[i]])$coefficients[2, 1]))
    }
}
all_transition_age_unc_cont=all_transition_age_unc_cont^2

if (number_of_models>=2){
  for (i in 1:number_of_models){
    for (j in extents[[i]]:extents[[i+1]]){
      if (j==extents[1]){
        CFCS_table[(j-extents[1]+1),6]=0
      }else{
        if (i==1){
          CFCS_table[(j-extents[1]+1),6]=model_transition_age_catch[(j-extents[1]+1)]*((-decay_const/summary(list_of_models[[i]])$coefficients[2, 1])*sqrt(((-summary(list_of_models[[i]])$coefficients[2, 2])/(summary(list_of_models[[i]])$coefficients[2,1]))^2+(decay_const_uncer/decay_const)^2)/(-decay_const/summary(list_of_models[[i]])$coefficients[2, 1]))
        }else{
          CFCS_table[(j-extents[1]+1),6]=sqrt(sum(all_transition_age_unc_cont[1:(i-1)])+(model_transition_age_catch[(j-extents[1]+1)]*((-decay_const/summary(list_of_models[[i]])$coefficients[2, 1])*sqrt(((-summary(list_of_models[[i]])$coefficients[2, 2])/(summary(list_of_models[[i]])$coefficients[2,1]))^2+(decay_const_uncer/decay_const)^2)/(-decay_const/summary(list_of_models[[i]])$coefficients[2, 1])))^2)
        }
      }
    }
  }
}else{ 
  for (i in extents[1]:extents[2]){
    if (i==extents[1]){
      CFCS_table[(i-extents[1]+1),6]=0
    }else{
      CFCS_table[(i-extents[1]+1),6]=CFCS_table[(i-extents[1]+1),4]*(MAR_final_err/MAR_final)
    }
  }
}

if(surface_active_zone_check==TRUE){
  has_known_age=readline(prompt = "Does one of your intervals have an independent age marker? (TRUE or FALSE) ")
    if (has_known_age==TRUE){
      known_age_depth=as.numeric(readline(prompt = "Enter marker depth (cm) "))
      known_age_CE=as.numeric(readline(prompt = "Enter marker age (CE) "))
      index <- which(CFCS_table$zi == known_age_depth)
      top_age=known_age_CE+CFCS_table[index,"Time yr"]
      CFCS_table[1,"Date (CE)"]=top_age
      for (i in 2:nrow(CFCS_table)){
        CFCS_table[i,"Date (CE)"]=top_age-CFCS_table[i,"Time yr"]
      }
    }
}


folder_name2 = "CFCS modeling data"
dir.create(folder_name2)
CFCS_table_path <- file.path(folder_name2, "CFCS_table.xlsx")
colnames(CFCS_table)[c(1, 2, 3, 5)] <- c("Mid Depth zi cm", "Mass Depth Section mi gcm^-2", "log-scale Excess Pb-210 ln(Ci) dpm/g", "Date CE")
write_xlsx(as.data.frame(CFCS_table), CFCS_table_path)
colnames(CFCS_table)[c(1, 2, 3, 5)] <- c("zi","mi","ln(Ci)", "Date (CE)")
# save the MAR/SAR rates and uncertainties to files
#MAR version
rownames_MAR_rates=c()
if (number_of_models==1){
  MAR_rates = array(data=NA, dim=c(1, 2),dimnames=list(NULL, c("MAR (g/(cm^2 yr)", "uncertainty")))
  MAR_rates[1,1]=-decay_const/summary(changepoint_detection_model)$coefficients[2, 1]
  MAR_rates[1,2]=MAR_rates[1,1]*sqrt((-summary(changepoint_detection_model)$coefficients[2, 2]/summary(changepoint_detection_model)$coefficients[2, 1])^2+(decay_const_uncer/decay_const)^2)
  rownames(MAR_rates)=c("single model")
}else{
  MAR_rates = array(data=NA, dim=c(length(list_of_models), 2),dimnames=list(NULL, c("MAR (g/(cm^2 yr)", "uncertainty")))
  for ( i in 1:length(list_of_models)){
    MAR_rates[i,1]=-decay_const/summary(list_of_models[[i]])$coefficients[2, 1]
    MAR_rates[i,2]=MAR_rates[i,1]*sqrt((-summary(list_of_models[[i]])$coefficients[2, 2]/summary(list_of_models[[i]])$coefficients[2, 1])^2+(decay_const_uncer/decay_const)^2)
    rownames_MAR_rates[i]=paste("Model",i)
  }
  rownames(MAR_rates)=rownames_MAR_rates
}
MAR_rates_path <- file.path(folder_name2, "MAR_rates.xlsx")
if (number_of_models>1){
  MAR_final_slope=summary(list_of_models[[length(list_of_models)]])$coefficients[2, "Estimate"]
  MAR_final_slope_err=summary(list_of_models[[length(list_of_models)]])$coefficients[2, 2]
}else{
  MAR_final_slope=summary(changepoint_detection_model)$coefficients[2, "Estimate"]
  MAR_final_slope_err=summary(changepoint_detection_model)$coefficients[2, 2]
}

write_xlsx(as.data.frame(MAR_rates), MAR_rates_path)
MAR_final=(-decay_const/MAR_final_slope)
MAR_final_err=(-decay_const/MAR_final_slope)*sqrt((MAR_final_slope_err/MAR_final_slope)^2+(decay_const_uncer/decay_const)^2)
#SAR version
rownames_SAR_rates=c()
if (number_of_models==1){
  SAR_rates = array(data=NA, dim=c(1, 2),dimnames=list(NULL, c("SAR (cm/yr)", "uncertainty")))
  SAR_rates[1,1]=-decay_const/summary(changepoint_detection_model_SAR)$coefficients[2, 1]
  SAR_rates[1,2]=SAR_rates[1,1]*sqrt((-summary(changepoint_detection_model_SAR)$coefficients[2, 2]/summary(changepoint_detection_model_SAR)$coefficients[2, 1])^2+(decay_const_uncer/decay_const)^2)
  rownames(SAR_rates)=c("single model")
}else{
  SAR_rates = array(data=NA, dim=c(length(list_of_models_SAR), 2),dimnames=list(NULL, c("SAR (cm/yr)", "uncertainty")))
  for ( i in 1:length(list_of_models)){
    SAR_rates[i,1]=-decay_const/summary(list_of_models_SAR[[i]])$coefficients[2, 1]
    SAR_rates[i,2]=SAR_rates[i,1]*sqrt((-summary(list_of_models_SAR[[i]])$coefficients[2, 2]/summary(list_of_models_SAR[[i]])$coefficients[2, 1])^2+(decay_const_uncer/decay_const)^2)
    rownames_SAR_rates[i]=paste("Model",i)
  }
  rownames(SAR_rates)=rownames_SAR_rates
}
SAR_rates_path <- file.path(folder_name2, "SAR_rates.xlsx")
write_xlsx(as.data.frame(SAR_rates), SAR_rates_path)

#ccSAR version
if (exists("zi_correction_table")){
  rownames_ccSAR_rates=c()
  if (number_of_models==1){
    ccSAR_rates = array(data=NA, dim=c(1, 2),dimnames=list(NULL, c("ccSAR (cm/yr)", "uncertainty")))
    ccSAR_rates[1,1]=-decay_const/summary(changepoint_detection_model_ccSAR)$coefficients[2, 1]
    ccSAR_rates[1,2]=ccSAR_rates[1,1]*sqrt((-summary(changepoint_detection_model_ccSAR)$coefficients[2, 2]/summary(changepoint_detection_model_ccSAR)$coefficients[2, 1])^2+(decay_const_uncer/decay_const)^2)
    rownames(ccSAR_rates)=c("single model")
  }else{
    ccSAR_rates = array(data=NA, dim=c(length(list_of_models_ccSAR), 2),dimnames=list(NULL, c("ccSAR (cm/yr)", "uncertainty")))
    for ( i in 1:length(list_of_models)){
      ccSAR_rates[i,1]=-decay_const/summary(list_of_models_ccSAR[[i]])$coefficients[2, 1]
      ccSAR_rates[i,2]=ccSAR_rates[i,1]*sqrt((-summary(list_of_models_ccSAR[[i]])$coefficients[2, 2]/summary(list_of_models_ccSAR[[i]])$coefficients[2, 1])^2+(decay_const_uncer/decay_const)^2)
      rownames_ccSAR_rates[i]=paste("Model",i)
    }
    rownames(ccSAR_rates)=rownames_ccSAR_rates
  }
  ccSAR_rates_path <- file.path(folder_name2, "ccSAR_rates.xlsx")
  write_xlsx(as.data.frame(ccSAR_rates), ccSAR_rates_path)
}

##save changepoint models, if they are used
if (number_of_models==1){
  MAR_model_summary <- data.frame(Model = character(),
                                  Coefficient = character(),
                                  Estimate = numeric(),
                                  Std.Error = numeric(),
                                  t.value = numeric(),
                                  Pr.t = numeric(),
                                  Adj.R.squared = numeric(),
                                  stringsAsFactors = FALSE)
  
  # Check if the model is valid
  if (!is.null(changepoint_detection_model)) {
    summary_model <- summary(changepoint_detection_model)
    
    # Extract coefficients
    coefficients <- summary_model$coefficients
    for (i in 1:nrow(coefficients)) {
      MAR_model_summary <- rbind(MAR_model_summary, 
                                 data.frame(Model = "Single Model",
                                            Coefficient = ifelse(rownames(coefficients)[i] == "(Intercept)", "Intercept", "Slope"),  # CHANGED
                                            Estimate = coefficients[i, "Estimate"],
                                            Std.Error = coefficients[i, "Std. Error"],
                                            t.value = coefficients[i, "t value"],
                                            Pr.t = coefficients[i, "Pr(>|t|)"],
                                            Adj.R.squared = summary_model$adj.r.squared,  # Add R-squared value
                                            stringsAsFactors = FALSE))
    }
    
    # If coefficients were extracted, add the model name
    if (nrow(MAR_model_summary) > 0) {
      MAR_model_summary$Model <- "Changepoint Detection Model"
    }
  } else {
    cat("The changepoint_detection_model is NULL or invalid.\n")
  }
  
  # Write the summary to an Excel file
  write_xlsx(MAR_model_summary, path = "CFCS modeling data/MAR_model_summary.xlsx")
  
  #now for the SAR version
  
  SAR_model_summary <- data.frame(Model = character(),
                                  Coefficient = character(),
                                  Estimate = numeric(),
                                  Std.Error = numeric(),
                                  t.value = numeric(),
                                  Pr.t = numeric(),
                                  Adj.R.squared = numeric(),
                                  stringsAsFactors = FALSE)
  
  # Check if the model is valid
  if (!is.null(changepoint_detection_model_SAR)) {
    summary_model <- summary(changepoint_detection_model_SAR)
    
    # Extract coefficients
    coefficients <- summary_model$coefficients
    for (i in 1:nrow(coefficients)) {
      SAR_model_summary <- rbind(SAR_model_summary, 
                                 data.frame(Model = "Single Model",
                                            Coefficient = ifelse(rownames(coefficients)[i] == "(Intercept)", "Intercept", "Slope"),  # CHANGED
                                            Estimate = coefficients[i, "Estimate"],
                                            Std.Error = coefficients[i, "Std. Error"],
                                            t.value = coefficients[i, "t value"],
                                            Pr.t = coefficients[i, "Pr(>|t|)"],
                                            Adj.R.squared = summary_model$adj.r.squared,  # Add R-squared value
                                            stringsAsFactors = FALSE))
    }
    
    # If coefficients were extracted, add the model name
    if (nrow(SAR_model_summary) > 0) {
      SAR_model_summary$Model <- "Changepoint Detection Model"
    }
  } else {
    cat("The changepoint_detection_model is NULL or invalid.\n")
  }
  
  # Write the summary to an Excel file
  write_xlsx(SAR_model_summary, path = "CFCS modeling data/SAR_model_summary.xlsx")

  #ccSAR version
  if (exists("zi_correction_table")){
    ccSAR_model_summary <- data.frame(Model = character(),
                                      Coefficient = character(),
                                      Estimate = numeric(),
                                      Std.Error = numeric(),
                                      t.value = numeric(),
                                      Pr.t = numeric(),
                                      Adj.R.squared = numeric(),
                                      stringsAsFactors = FALSE)
    
    # Check if the model is valid
    if (!is.null(changepoint_detection_model_ccSAR)) {
      ccsummary_model <- summary(changepoint_detection_model_ccSAR)
      
      # Extract coefficients
      coefficients <- ccsummary_model$coefficients
      for (i in 1:nrow(coefficients)) {
        ccSAR_model_summary <- rbind(ccSAR_model_summary, 
                                     data.frame(Model = "Single Model",
                                                Coefficient = rownames(coefficients)[i],
                                                Estimate = coefficients[i, "Estimate"],
                                                Std.Error = coefficients[i, "Std. Error"],
                                                t.value = coefficients[i, "t value"],
                                                Pr.t = coefficients[i, "Pr(>|t|)"],
                                                Adj.R.squared = summary_model$adj.r.squared,  # Add R-squared value
                                                stringsAsFactors = FALSE))
      }
      
      # If coefficients were extracted, add the model name
      if (nrow(ccSAR_model_summary) > 0) {
        ccSAR_model_summary$Model <- "Changepoint Detection Model"
      }
    } else {
      cat("The changepoint_detection_model is NULL or invalid.\n")
    }
    
    # Write the summary to an Excel file
    write_xlsx(ccSAR_model_summary, path = "CFCS modeling data/ccSAR_model_summary.xlsx")
  }
  
  
  }

## MAR lm summary save
if(number_of_models>1){
  MAR_model_summaries <- data.frame(Model = character(),
                                    Coefficient = character(),
                                    Estimate = numeric(),
                                    Std.Error = numeric(),
                                    t.value = numeric(),
                                    Pr.t = numeric(),
                                    Adj.R.squared = numeric(),
                                    stringsAsFactors = FALSE)
  
  # Extract information from each model and store it in the data frame
  for (model_index in seq_along(list_of_models)) {
    model <- list_of_models[[model_index]]
    
    # Check if the model is valid
    if (!is.null(model)) {
      summary_model <- summary(model)
      
      # Extract coefficients
      coefficients <- summary_model$coefficients
      for (i in 1:nrow(coefficients)) {
        MAR_model_summaries <- rbind(MAR_model_summaries, 
                                     data.frame(Model = paste("Model", model_index, "Indicies",extents[model_index] , "-", extents[model_index+1]),  # Adjusting index for naming
                                                Coefficient = ifelse(rownames(coefficients)[i] == "(Intercept)", "Intercept", "Slope"),  # CHANGED
                                                Estimate = coefficients[i, "Estimate"],
                                                Std.Error = coefficients[i, "Std. Error"],
                                                t.value = coefficients[i, "t value"],
                                                Pr.t = coefficients[i, "Pr(>|t|)"],
                                                Adj.R.squared = summary_model$adj.r.squared,  # Add R-squared value
                                                stringsAsFactors = FALSE))
      }
    } else {
      cat("Model at index", model_index, "is NULL or invalid.\n")
    }
  }
  write_xlsx(MAR_model_summaries, path = "CFCS modeling data/MAR_model_summaries.xlsx")
}
## SAR model summaries
## SAR lm summary save
if(number_of_models>1){
  SAR_model_summaries <- data.frame(Model = character(),
                                    Coefficient = character(),
                                    Estimate = numeric(),
                                    Std.Error = numeric(),
                                    t.value = numeric(),
                                    Pr.t = numeric(),
                                    Adj.R.squared = numeric(),
                                    stringsAsFactors = FALSE)
  
  # Extract information from each model and store it in the data frame
  for (model_index in seq_along(list_of_models_SAR)) {
    model <- list_of_models_SAR[[model_index]]
    
    # Check if the model is valid
    if (!is.null(model)) {
      summary_model <- summary(model)
      
      # Extract coefficients
      coefficients <- summary_model$coefficients
      for (i in 1:nrow(coefficients)) {
        SAR_model_summaries <- rbind(SAR_model_summaries, 
                                     data.frame(Model = paste("Model", model_index, "Indicies",extents[model_index] , "-", extents[model_index+1]),  # Adjusting index for naming
                                                Coefficient = ifelse(rownames(coefficients)[i] == "(Intercept)", "Intercept", "Slope"),  # CHANGED
                                                Estimate = coefficients[i, "Estimate"],
                                                Std.Error = coefficients[i, "Std. Error"],
                                                t.value = coefficients[i, "t value"],
                                                Pr.t = coefficients[i, "Pr(>|t|)"],
                                                Adj.R.squared = summary_model$adj.r.squared,  # Add R-squared value
                                                stringsAsFactors = FALSE))
      }
    } else {
      cat("Model at index", model_index, "is NULL or invalid.\n")
    }
  }
  write_xlsx(SAR_model_summaries, path = "CFCS modeling data/SAR_model_summaries.xlsx")
}

#ccSAR version
if (exists("zi_correction_table")){
  ## ccSAR model summaries
  ## ccSAR lm summary save
  if(number_of_models>1){
    ccSAR_model_summaries <- data.frame(Model = character(),
                                        Coefficient = character(),
                                        Estimate = numeric(),
                                        Std.Error = numeric(),
                                        t.value = numeric(),
                                        Pr.t = numeric(),
                                        Adj.R.squared = numeric(),
                                        stringsAsFactors = FALSE)
    
    # Extract information from each model and store it in the data frame
    for (model_index in seq_along(list_of_models_ccSAR)) {
      model <- list_of_models_ccSAR[[model_index]]
      
      # Check if the model is valid
      if (!is.null(model)) {
        summary_model <- summary(model)
        
        # Extract coefficients
        coefficients <- summary_model$coefficients
        for (i in 1:nrow(coefficients)) {
          ccSAR_model_summaries <- rbind(ccSAR_model_summaries, 
                                         data.frame(Model = paste("Model", model_index, "Indicies",extents[model_index] , "-", extents[model_index+1]),  # Adjusting index for naming
                                                    Coefficient = ifelse(rownames(coefficients)[i] == "(Intercept)", "Intercept", "Slope"),  # CHANGED
                                                    Estimate = coefficients[i, "Estimate"],
                                                    Std.Error = coefficients[i, "Std. Error"],
                                                    t.value = coefficients[i, "t value"],
                                                    Pr.t = coefficients[i, "Pr(>|t|)"],
                                                    Adj.R.squared = summary_model$adj.r.squared,  # Add R-squared value
                                                    stringsAsFactors = FALSE))
        }
      } else {
        cat("Model at index", model_index, "is NULL or invalid.\n")
      }
    }
    write_xlsx(ccSAR_model_summaries, path = "CFCS modeling data/ccSAR_model_summaries.xlsx")
  }
}

# plot saves
# SAR plot save
pdf(
  file.path(folder_name2, "SAR_models.pdf"),
  width = 8,
  height = 6
)

if (number_of_models > 1) {
  
  SAR_models_plot <- ggplot(
    CFCS_table,
    aes(x = zi, y = `ln(Ci)`)
  ) +
    geom_point(
      color = "blue"
    ) +
    labs(
      title = "SAR Modeling",
      x = "Mid Depth (cm)",
      y = "log-scale Pb-210 Activity (dpm/g)"
    ) +
    geom_line(
      data = lines_df_SAR,
      aes(
        x = x,
        y = y,
        group = group,
        color = group
      ),
      linewidth = 1.2,
      show.legend = TRUE
    ) +
    labs(
      color = "Model"
    ) +
    dynamic_scale(
      CFCS_table,
      "zi"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  
  print(SAR_models_plot)
  
} else {
  
  print(
    SAR_model_plot +
      geom_line(
        data = single_line_df_SAR,
        aes(x = zi, y = fitted_y),
        inherit.aes = FALSE,
        color = "red",
        linewidth = 1
      )
  )
}

dev.off()
#ccSAR version
if (exists("zi_correction_table")){
  # ccSAR plot save
  pdf(
    file.path(folder_name2, "ccSAR_models.pdf"),
    width = 8,
    height = 6
  )
  
  if (number_of_models > 1) {
    
    ccSAR_models_plot <- ggplot(
      CFCS_table,
      aes(x = CFCS_table[, "Uncompacted Mid-Depth (cm)"], y = `ln(Ci)`)
    ) +
      geom_point(
        color = "blue"
      ) +
      labs(
        title = "Uncompacted SAR Modeling",
        x = "Uncompacted Mid Depth (cm)",
        y = "log-scale Pb-210 Activity (dpm/g)"
      ) +
      geom_line(
        data = lines_df_ccSAR,
        aes(
          x = x,
          y = y,
          group = group,
          color = group
        ),
        linewidth = 1.2,
        show.legend = TRUE
      ) +
      labs(
        color = "Model"
      ) +
      dynamic_scale(
        CFCS_table,
        "Uncompacted Mid-Depth (cm)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5)
      )
    
    print(ccSAR_models_plot)
    
  } else {
    
    print(
      ccSAR_model_plot +
        geom_line(
          data = single_line_df_ccSAR,
          aes(x = x, y = fitted_y),
          inherit.aes = FALSE,
          color = "red",
          linewidth = 1
        )
    )
  }
  
  dev.off()
}

# MAR plot save
pdf(
  file.path(folder_name2, "MAR_models.pdf"),
  width = 8,
  height = 6
)

if (number_of_models > 1) {
  
  MAR_models_plot <- ggplot(
    CFCS_table,
    aes(x = mi, y = `ln(Ci)`)
  ) +
    geom_point(
      color = "blue"
    ) +
    labs(
      title = "MAR Modeling",
      x = "Cumulative Mass Depth (g/cm^2)",
      y = "log-scale Pb-210 Activity (dpm/g)"
    ) +
    geom_line(
      data = lines_df,
      aes(
        x = x,
        y = y,
        group = group,
        color = group
      ),
      linewidth = 1.2,
      show.legend = TRUE
    ) +
    labs(
      color = "Model"
    ) +
    dynamic_scale(
      CFCS_table,
      "mi"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  
  print(MAR_models_plot)
  
} else {
  
  print(
    MAR_model_plot +
      geom_line(
        data = single_line_df,
        aes(x = mi, y = fitted_y),
        inherit.aes = FALSE,
        color = "red",
        linewidth = 1
      )
  )
}

dev.off()

# Age-depth
pdf(
  file.path(folder_name2, "age_depth.pdf"),
  width = 8,
  height = 6
)

age_depth_plot <- ggplot(
  CFCS_table,
  aes(x = `Date (CE)`, y = `zi`)
) +
  geom_point(
    color = "blue"
  ) +
  geom_errorbar(
    aes(
      xmin = `Date (CE)` - `Time Uncertainty (yr)`,
      xmax = `Date (CE)` + `Time Uncertainty (yr)`,
      y = `zi`
    ),
    width = 0.2,
    color = "black"
  ) +
  labs(
    title = "Age vs Depth Model",
    y = "Mid Depth (cm)",
    x = "Year (CE)"
  ) +
  dynamic_scale(
    CFCS_table,
    "Date (CE)",
    reverse = TRUE
  ) +
  dynamic_scale(
    CFCS_table,
    "zi",
    axis = "y",
    reverse = TRUE
  ) +
  theme_minimal()

print(age_depth_plot)

dev.off()

#ccSAR version
if (exists("zi_correction_table")){
  # Age-uncompacted-depth
  pdf(
    file.path(folder_name2, "age_unccompacted_depth.pdf"),
    width = 8,
    height = 6
  )
  
  age_depth_plot <- ggplot(
    CFCS_table,
    aes(x = `Date (CE)`, y = `Uncompacted Mid-Depth (cm)`)
  ) +
    geom_point(
      color = "blue"
    ) +
    geom_errorbar(
      aes(
        xmin = `Date (CE)` - `Time Uncertainty (yr)`,
        xmax = `Date (CE)` + `Time Uncertainty (yr)`,
        y = `zi`
      ),
      width = 0.2,
      color = "black"
    ) +
    labs(
      title = "Age vs Depth Model",
      y = "Uncompacted Mid Depth (cm)",
      x = "Year (CE)"
    ) +
    dynamic_scale(
      CFCS_table,
      "Date (CE)",
      reverse = TRUE
    ) +
    dynamic_scale(
      CFCS_table,
      "Uncompacted Mid-Depth (cm)",
      axis = "y",
      reverse = TRUE
    ) +
    theme_minimal()
  
  print(age_depth_plot)
  
  dev.off()
}

# Age-mass depth
pdf(
  file.path(folder_name2, "age_mass_depth.pdf"),
  width = 8,
  height = 6
)

age_mass_depth_plot <- ggplot(
  CFCS_table,
  aes(x = `Date (CE)`, y = `mi`)
) +
  geom_point(
    color = "blue"
  ) +
  geom_errorbar(
    aes(
      xmin = `Date (CE)` - `Time Uncertainty (yr)`,
      xmax = `Date (CE)` + `Time Uncertainty (yr)`,
      y = `mi`
    ),
    width = 0.2,
    color = "black"
  ) +
  labs(
    title = "Age vs Mass Depth Model",
    y = "Cumulative Mass Depth (g/cm^2)",
    x = "Year (CE)"
  ) +
  dynamic_scale(
    CFCS_table,
    "Date (CE)",
    reverse = TRUE
  ) +
  dynamic_scale(
    CFCS_table,
    "mi",
    axis = "y",
    reverse = TRUE
  ) +
  theme_minimal()

print(age_mass_depth_plot)

dev.off()

if (surface_active_zone_check==TRUE){
  print("CFCS model assumptions violated! This core has a surface active zone, only the relative ages of the sediment layers below the surface active layer could be calculated, not absolute ages")
  CFCS_table[,5]=NA
}

if (exists("hit_background") && !hit_background) {
  print("CFCS model assumptions violated! Supported Pb-210 activity is unknown, CFCS is not a valid model")
}
inventory_completed=FALSE
missing_inventory_method=readline(prompt = "Can CFCS be used to complete missing inventory? (TRUE or FALSE) ")
if (missing_inventory_method==TRUE){
  missing_inventory=(concentrations_table[(nrow(concentrations_table)),6]*(-decay_const/MAR_final_slope))/decay_const
  missing_inventory_uncer=missing_inventory*sqrt((decay_const_uncer/decay_const)^2+(MAR_final_err/MAR_final)^2+(concentrations_table[nrow(concentrations_table),7]/concentrations_table[nrow(concentrations_table),6])^2)
inventory_completed=TRUE
  }else{
ref_date_check=readline(prompt = "Do you have an independent chronometer reference date for a sediment layer? (TRUE or FALSE) ")  
if (ref_date_check==TRUE){
  ref_date=as.numeric(readline(prompt = "Enter reference age (CE) "))
  ref_depth=as.numeric(readline(prompt = "Enter reference date depth (cm) "))
  # Locate the independently dated reference layer
  Aj_idx <- which(!is.na(concentrations_table[, 1]) & concentrations_table[, 1] == ref_depth)
  # Reference-layer age
  delta_t <- numeric_date - ref_date
  D <- exp(delta_t * decay_const) - 1
  # Identify rows above and below the reference layer.
  upper_idx <- if (Aj_idx > 1) {
    1:(Aj_idx - 1)
  } else {
    integer(0)
  }
  
  lower_idx <- if (Aj_idx < nrow(concentrations_table)) {
    (Aj_idx + 1):nrow(concentrations_table)
  } else {
    integer(0)
  }
  # Reference-layer inventory and uncertainty
  reference_inventory <- concentrations_table[Aj_idx, 10]
  reference_error <- concentrations_table[Aj_idx, 11]
  # Measured inventory above the reference layer: delta A
  A_delta <- sum(concentrations_table[upper_idx, 10],na.rm = TRUE)+0.5*reference_inventory
  # Uncertainty in A_delta
  A_delta_error <- sqrt(sum(concentrations_table[upper_idx, 11]^2,na.rm = TRUE)+(0.5*reference_error)^2)
  # Measured inventory below the reference layer
  Aj_measured_fraction <- sum(concentrations_table[lower_idx, 10],na.rm = TRUE)+0.5*reference_inventory
  # Uncertainty in measured inventory below reference
  Aj_measured_error <- sqrt(sum(concentrations_table[lower_idx, 11]^2,na.rm = TRUE)+(0.5 * reference_error)^2)
  # Calculate complete inventory below the reference horizon
  Aj_calculated <- A_delta / D
  # Uncertainty in Aj_calculated
  Aj_calculated_error <-A_delta_error / D
  # Covariance between Aj_calculated and Aj_measured_fraction
  covariance_calc_measured <-reference_error^2 / (4 * D)
  # Missing inventory
  missing_inventory <-Aj_calculated - Aj_measured_fraction
  missing_inventory_variance <-Aj_calculated_error^2 +Aj_measured_error^2 -2 * covariance_calc_measured
  # Prevent tiny negative values from floating-point round-off
  missing_inventory_variance <-max(0, missing_inventory_variance)
  missing_inventory_uncer <-sqrt(missing_inventory_variance)
  
inventory_completed=TRUE
  }else{
ind_MAR_check=readline(prompt = "Do you have an independent estimate of MAR at the bottom of the sediment core? (TRUE or FALSE) ")
if (ind_MAR_check==TRUE){
  independent_MAR_num=as.numeric(readline(prompt = "Enter MAR (g/(cm^2 yr)) "))
  independent_MAR_num_uncer=as.numeric(readline(prompt = "Enter MAR uncertainty (g/(cm^2 yr)) "))
  missing_inventory=(concentrations_table[(nrow(concentrations_table)),6]*(independent_MAR_num))/decay_const
  missing_inventory_uncer=missing_inventory*sqrt((decay_const_uncer/decay_const)^2+(independent_MAR_num_uncer/independent_MAR_num)^2+(concentrations_table[nrow(concentrations_table),7]/concentrations_table[nrow(concentrations_table),6])^2)
inventory_completed=TRUE
  }
}
}
if (inventory_completed==FALSE){
  missing_inventory=0
  missing_inventory_uncer=0
}

finish_CFCS=readline(prompt = "CFCS modeling complete. Continue to CF? (True or False) ")
#start of CF section
if (finish_CFCS==TRUE){
#create CF table and begin populating it
CF_table=array(data=NA, dim=c(nrow(mass_table),17), dimnames=list(NULL,c("Top of Interval z(i) cm", "Mid Depth zi cm","Inventory delta Ai dpm/cm^2","Inventory Uncertainty u(delta Ai) dpm/cm^2","Inventory Below A(i) dpm/cm^2","Inventory Below Uncertainty u(A(i)) dpm/cm^2","Age t(i) yr","Age Uncertainty u(t(i)) yr","Year CE","Mass Accumulation Rate r(i) g/(cm^2 yr)","Mass Accumulation Rate Uncertainty u(r(i)) g/(cm^2 yr)","Dry Bulk Density Section ri g/cm^3","Dry Bulk Density Section Uncertainty u(ri) g/cm^3","Dry Bulk Density Layer r(i) g/cm^3","Dry Bulk Density Layer Uncertainty u(r(i)) g/cm^3","Sediment Accumulation Rate s(i) cm/yr","Sediment Accumulation Rate Uncertainty U(s(i)) cm/yr")))
colnames(CF_table) <- gsub("delta", "\u0394", colnames(CF_table))
CF_table[,1:2]=mass_table[,1:2]
CF_table[1:nrow(concentrations_table),3]=concentrations_table[,10]
CF_table[1:nrow(concentrations_table),4]=concentrations_table[,11]
CF_table[,12]=mass_table[,4]
CF_table[,13]=mass_table[,5]
#inventory below and uncertainty
for (i in (nrow(CF_table)):1){
  if (i %% 2 ==1){
    if (i==nrow(CF_table)){
      CF_table[nrow(CF_table),5]=missing_inventory
      CF_table[nrow(CF_table),6]=missing_inventory_uncer
    }else{
      CF_table[i,5]=CF_table[i+2,5]+CF_table[i+1,3]
      CF_table[i, 6] = sqrt(sum(CF_table[i + 2, 6]^2, CF_table[i+1, 4]^2))
    }
  }
}
#ti, year, MAR
for (i in 1:(nrow(CF_table)-2)){
  if (i %% 2 ==1){
    CF_table[i,7]=log(CF_table[1,5]/CF_table[i,5])/decay_const
    CF_table[i,9]=numeric_date-CF_table[i,7]
    CF_table[i,10]=(decay_const*CF_table[i,5])/concentrations_table[i,8]
  }
}
# Age uncertainty loop
for (i in 1:(nrow(CF_table) - 2)) {
  if (i %% 2 == 1) {
    # Calculate the components of the formula
    term1 <- (CF_table[i, 7] * decay_const_uncer)^2
    term2 <- (CF_table[1, 6] / CF_table[1, 5])^2  # Use the first row values
    term3 <- (1 - 2 * (CF_table[i, 5] / CF_table[1, 5])) * (CF_table[i, 6] / CF_table[i, 5])^2
    
    # Combine the terms and calculate the final value
    CF_table[i, 8] <- sqrt(term1 + term2 + term3) / decay_const
  }
}
# MAR uncertainty
for (i in 1:(nrow(CF_table)-2)){
  if (i %% 2 ==1){
    CF_table[i,11]=CF_table[i,10]*sqrt((CF_table[i,6]/CF_table[i,5])^2+(1-((CF_table[i,5]-CF_table[(i+2),5])/CF_table[i,5]))*(concentrations_table[i,9]/concentrations_table[i,8])^2)
  }
}
# surface DBD determination
has_DBD_0=readline(prompt="Do you have a measured DBD surface value? (TRUE or FALSE)")

if (has_DBD_0==TRUE){
  DBD_layer=as.numeric(readline(prompt="Enter DBD (g/cm^3)"))
  DBD_layer_uncer=as.numeric(readline(prompt="Enter DBD Uncertainty (g/cm^3)"))
}else{

DBD_layer_table=array(data=NA, dim=c(10,2), dimnames=list(NULL, c("zi","DBDi")))

# Construct the first 10 interval records
DBD_layer_table <- data.frame(
  zi = mass_table[seq(2, 20, by = 2), 2],
  DBDi = NA_real_
)

DBD <- mass_table[
  seq(2, 20, by = 2),
  "Dry Bulk Density DBD gcm^-3"
]

# Only take logarithms of finite, positive activities
valid_DBD <- is.finite(DBD) & DBD > 0

DBD_layer_table$DBDi[valid_DBD] <- DBD[valid_DBD]

# Preserve the original-measurement rule
for (j in seq_len(nrow(DBD_layer_table))) {
  if (is.na(original_section_DBD_evens[j]) ||
      !original_section_DBD_evens[j]) {
    DBD_layer_table$DBDi[j] <- NA_real_
  }
}

DBD_layer_model_list <- list()
model_summaries_surf_DBD <- list()

for (extent in 3:10) {
  
  candidate <- DBD_layer_table[seq_len(extent), ]
  
  # Identify valid observations
  valid <- is.finite(candidate[, "zi"]) &
    is.finite(candidate[, "DBDi"])
  
  fit_data <- candidate[valid, , drop = FALSE]
  
  # Use the same columns for every candidate model
  result <- data.frame(
    Model = paste0("Indices 1-", extent),
    N = nrow(fit_data),
    Used_indices = paste(which(valid), collapse = ", "),
    Intercept = NA_real_,
    Slope = NA_real_,
    Intercept_SE = NA_real_,
    Slope_SE = NA_real_,
    R.squared = NA_real_,
    Adj.R.squared = NA_real_,
    Status = "Insufficient valid observations",
    stringsAsFactors = FALSE
  )
  
  # Require at least 3 observations and distinct predictor values
  if (nrow(fit_data) >= 3 &&
      length(unique(fit_data$zi)) >= 2) {
    
    model <- lm(DBDi ~ zi, data = fit_data)
    
    DBD_layer_model_list[[as.character(extent - 2)]] <- model
    
    sm <- summary(model)
    coefs <- sm$coefficients
    
    result$Intercept <- unname(coef(model)[1])
    result$Slope <- unname(coef(model)[2])
    result$Intercept_SE <- coefs[1, "Std. Error"]
    result$Slope_SE <- coefs[2, "Std. Error"]
    result$R.squared <- sm$r.squared
    result$Adj.R.squared <- sm$adj.r.squared
    result$Status <- "Fitted"
  }
  
  model_summaries_surf_DBD[[as.character(extent)]] <- result
}

# Combine all candidate summaries
model_summary_table <- do.call(
  rbind,
  model_summaries_surf_DBD
)
model_summaries_surf_DBD=model_summary_table
rownames(model_summaries_surf_DBD) <- NULL
print("Please determine surface layer DBD")
print(DBD_layer_table)
print(model_summaries_surf_DBD[,-(2:3)])
selection_complete=FALSE
while (selection_complete==FALSE){
  best_model_index=readline(prompt = "Select model (1-8) ")
  
  # Step 3: Retrieve the corresponding model
  best_model <- DBD_layer_model_list[[best_model_index]]
  
  DBD_layer <- summary(best_model)$coefficients[1,1]
  DBD_layer_uncer <- summary(best_model)$coefficients["(Intercept)", "Std. Error"]
  print("Surface DBD (dpm/g)")
  print(unname(DBD_layer))
  selection_complete=readline(prompt = "Accept model selection? (TRUE or FALSE) ")
}
}
#DBD loop
for (i in 1:(nrow(CF_table)-2)){
  if (i %% 2 ==1){
    if (i==1){
      CF_table[1,14]=DBD_layer
    }else{
      fraction =
        (CF_table[i,"Top of Interval z(i) cm"] -
           CF_table[(i-1),"Mid Depth zi cm"]) /
        (CF_table[(i+1),"Mid Depth zi cm"] -
           CF_table[(i-1),"Mid Depth zi cm"])
      
      CF_table[i,"Dry Bulk Density Layer r(i) g/cm^3"] =
        CF_table[(i-1),"Dry Bulk Density Section ri g/cm^3"] +
        ((CF_table[(i+1),"Dry Bulk Density Section ri g/cm^3"] -
            CF_table[(i-1),"Dry Bulk Density Section ri g/cm^3"])) * fraction
      
      CF_table[i,"Dry Bulk Density Layer Uncertainty u(r(i)) g/cm^3"] =
        sqrt((1-fraction)^2 * CF_table[(i-1),"Dry Bulk Density Section Uncertainty u(ri) g/cm^3"]^2 +
               (fraction)^2 * CF_table[(i+1),"Dry Bulk Density Section Uncertainty u(ri) g/cm^3"]^2)
    }
  }
}
# SAR loop
for (i in 1:(nrow(CF_table)-2)){
  if (i %% 2 ==1){
    CF_table[i,16]=CF_table[i,10]/CF_table[i,14]
  }
}
# DBD section uncertainty
for (i in 1:(nrow(CF_table)-2)){
  if (i %% 2 ==1){
    if (i==1){
      CF_table[i,15]=DBD_layer_uncer
    }else{
      CF_table[i,15]=sqrt((CF_table[(i+1),13])^2+(CF_table[(i-1),13])^2)/2
    }
  }
}

#SAR uncertainty loop
for (i in 1:(nrow(CF_table)-2)){
  if (i %% 2 ==1){
    CF_table[i,17]=CF_table[i,16]*sqrt((CF_table[i,11]/CF_table[i,10])^2+(CF_table[i,15]/CF_table[i,14])^2)
  }
}
if (exists("zi_correction_table")){
# ============================================================
# CF/CRS LAYER-LEVEL COMPACTION CORRECTION
# ============================================================
#
# Reconstructs depositional layer DBD and calculates
# compaction-corrected CF/CRS SAR.
#
# Original CF_table and mass_table remain unchanged.
#
# ============================================================


# ------------------------------------------------------------
# 1. VALIDATE REQUIRED INPUTS
# ------------------------------------------------------------

n_sections <- length(section_rows)

# CF/CRS infinitesimal layer rows
layer_rows <- seq(1, nrow(CF_table) - 2, by = 2)

n_layers <- length(layer_rows)

# Section midpoint depths
section_depths <- as.numeric(
  mass_table[section_rows, "Mid Depth zi cm"]
)

# Layer boundary depths
layer_depths <- as.numeric(
  CF_table[layer_rows, "Top of Interval z(i) cm"]
)

# Section porosity
phi_sections <- as.numeric(
  zi_correction_table[
    section_rows, "Porosity"
  ]
)

# ------------------------------------------------------------
# 2. EXTRACT ORIGINAL CF/CRS LAYER VALUES
# ------------------------------------------------------------

MAR_original <- as.numeric(
  CF_table[layer_rows, 10]
)

uMAR_original <- as.numeric(
  CF_table[layer_rows, 11]
)

DBD_layer_original <- as.numeric(
  CF_table[layer_rows, 14]
)

uDBD_layer_original <- as.numeric(
  CF_table[layer_rows, 15]
)

SAR_original <- as.numeric(
  CF_table[layer_rows, 16]
)

uSAR_original <- as.numeric(
  CF_table[layer_rows, 17]
)

Age_original <- as.numeric(
  CF_table[layer_rows, 7]
)

Year_original <- as.numeric(
  CF_table[layer_rows, 9]
)

if (any(!is.finite(MAR_original)) ||
    any(!is.finite(uMAR_original)) ||
    any(!is.finite(DBD_layer_original)) ||
    any(DBD_layer_original <= 0)) {
  stop("Invalid original CF/CRS MAR or layer DBD.")
}


# ------------------------------------------------------------
# 3. INTERPOLATE PRESENT-DAY LAYER POROSITY
# ------------------------------------------------------------
#
# Interpolation uses actual linear depth.
#
# At the surface:
# phi_layer(0) = phi_0
#
# At all other layer depths:
# Linear interpolation between section midpoints.
#
# W_phi maps:
# [phi_section(1), ..., phi_section(n), phi_0]
# onto the infinitesimal layer porosities.
#
# ------------------------------------------------------------

W_phi_layer <- matrix(
  0,
  nrow = n_layers,
  ncol = n_sections + 1
)

for (i in seq_len(n_layers)) {
  
  z <- layer_depths[i]
  
  # Surface layer
  
  if (abs(z) < 1e-10) {
    
    W_phi_layer[i, n_sections + 1] <- 1
    
  } else {
    
    # Find the two section midpoints bounding this depth
    
    j <- findInterval(z, section_depths)
    
    if (j < 1 || j >= n_sections) {
      stop(paste(
        "Layer depth falls outside the section midpoint range:",
        z
      ))
    }
    
    k <- j + 1
    
    fraction <- (
      z - section_depths[j]
    ) / (
      section_depths[k] - section_depths[j]
    )
    
    if (!is.finite(fraction) ||
        fraction < 0 || fraction > 1) {
      stop(paste(
        "Invalid porosity interpolation fraction at layer", i
      ))
    }
    
    W_phi_layer[i, j] <- 1 - fraction
    W_phi_layer[i, k] <- fraction
  }
}


# Calculate interpolated layer porosity

phi_layer <- as.numeric(
  W_phi_layer %*%
    c(phi_sections, phi_0)
)


# ------------------------------------------------------------
# 4. PROPAGATE LAYER POROSITY COVARIANCE
# ------------------------------------------------------------

Sigma_phi_layer <- W_phi_layer %*%
  Sigma_input %*%
  t(W_phi_layer)

Sigma_phi_layer <- (
  Sigma_phi_layer + t(Sigma_phi_layer)
) / 2

u_phi_layer <- sqrt(
  pmax(diag(Sigma_phi_layer), 0)
)


# Preserve covariance between layer porosity and phi_0

W_phi_augmented <- rbind(
  W_phi_layer,
  c(rep(0, n_sections), 1)
)

Sigma_phi_augmented <- W_phi_augmented %*%
  Sigma_input %*%
  t(W_phi_augmented)

Sigma_phi_augmented <- (
  Sigma_phi_augmented +
    t(Sigma_phi_augmented)
) / 2


# ------------------------------------------------------------
# 5. CONSTRUCT LAYER DBD COVARIANCE
# ------------------------------------------------------------
#
# The first layer uses the independently estimated
# present-day surface DBD.
#
# Remaining layer DBD values are interpolated between
# section midpoint DBD values.
#
# Where available, use the full section DBD covariance
# from the preceding DBD correction block.
#
# Otherwise, use the section DBD uncertainties as
# independent measurements.
#
# ------------------------------------------------------------

DBD_sections <- as.numeric(
  mass_table[section_rows,
             "Dry Bulk Density DBD gcm^-3"]
)

uDBD_sections <- as.numeric(
  mass_table[section_rows,
             "Dry Bulk Density Uncertainty u(DBD) gcm^-3"]
)


# Section DBD covariance

if (exists("Sigma_DBD_present")) {
  
  Sigma_DBD_sections <- Sigma_DBD_present
  
} else {
  
  Sigma_DBD_sections <- diag(
    uDBD_sections^2,
    nrow = n_sections,
    ncol = n_sections
  )
}


# Surface DBD and uncertainty

DBD_surface <- as.numeric(DBD_layer_original[1])

uDBD_surface <- as.numeric(uDBD_layer_original[1])

if (!is.finite(DBD_surface) ||
    !is.finite(uDBD_surface) ||
    DBD_surface <= 0 ||
    uDBD_surface < 0) {
  stop("Invalid surface DBD or uncertainty.")
}


# Layer DBD interpolation matrix
#
# Input order:
# [surface DBD, section DBD(1), ..., section DBD(n)]

W_DBD_layer <- matrix(
  0,
  nrow = n_layers,
  ncol = n_sections + 1
)

for (i in seq_len(n_layers)) {
  
  z <- layer_depths[i]
  
  if (abs(z) < 1e-10) {
    
    W_DBD_layer[i, 1] <- 1
    
  } else {
    
    j <- findInterval(z, section_depths)
    
    if (j < 1 || j >= n_sections) {
      stop(paste(
        "Layer depth falls outside DBD interpolation range:",
        z
      ))
    }
    
    k <- j + 1
    
    fraction <- (
      z - section_depths[j]
    ) / (
      section_depths[k] - section_depths[j]
    )
    
    W_DBD_layer[i, j + 1] <- 1 - fraction
    W_DBD_layer[i, k + 1] <- fraction
  }
}


# Input DBD covariance

Sigma_DBD_layer_input <- matrix(
  0,
  nrow = n_sections + 1,
  ncol = n_sections + 1
)

Sigma_DBD_layer_input[1, 1] <- uDBD_surface^2

Sigma_DBD_layer_input[
  2:(n_sections + 1),
  2:(n_sections + 1)
] <- Sigma_DBD_sections


# Layer DBD covariance

Sigma_DBD_layer_present <- W_DBD_layer %*%
  Sigma_DBD_layer_input %*%
  t(W_DBD_layer)

Sigma_DBD_layer_present <- (
  Sigma_DBD_layer_present +
    t(Sigma_DBD_layer_present)
) / 2

uDBD_layer_present <- sqrt(
  pmax(diag(Sigma_DBD_layer_present), 0)
)


# Check interpolation against existing CF layer DBD

DBD_layer_check <- as.numeric(
  W_DBD_layer %*%
    c(DBD_surface, DBD_sections)
)

if (!isTRUE(all.equal(
  DBD_layer_check,
  DBD_layer_original,
  tolerance = 1e-8
))) {
  stop(
    "Reconstructed layer DBD does not match existing CF_table."
  )
}


# ------------------------------------------------------------
# 6. CALCULATE DEPOSITIONAL LAYER DBD
# ------------------------------------------------------------

DBD_layer_depositional <- DBD_layer_original *
  (1 - phi_0) / (1 - phi_layer)


# ------------------------------------------------------------
# 7. CONSTRUCT COMBINED INPUT COVARIANCE
# ------------------------------------------------------------
#
# Input vector:
#
# [DBD_layer(1:n),
#  phi_layer(1:n),
#  phi_0]
#
# DBD and porosity measurement errors are assumed
# independent.
#
# ------------------------------------------------------------

Sigma_layer_input <- matrix(
  0,
  nrow = 2 * n_layers + 1,
  ncol = 2 * n_layers + 1
)

DBD_indices <- seq_len(n_layers)

phi_indices <- (
  n_layers + 1
):(2 * n_layers + 1)

Sigma_layer_input[
  DBD_indices,
  DBD_indices
] <- Sigma_DBD_layer_present

Sigma_layer_input[
  phi_indices,
  phi_indices
] <- Sigma_phi_augmented

Sigma_layer_input <- (
  Sigma_layer_input + t(Sigma_layer_input)
) / 2


# ------------------------------------------------------------
# 8. JACOBIAN FOR DEPOSITIONAL LAYER DBD
# ------------------------------------------------------------
#
# DBD* = DBD * (1-phi_0)/(1-phi_layer)
#
# Derivatives:
#
# dDBD*/dDBD = (1-phi_0)/(1-phi_layer)
#
# dDBD*/dphi_layer = DBD*/(1-phi_layer)
#
# dDBD*/dphi_0 = -DBD*/(1-phi_0)
#
# ------------------------------------------------------------

J_layer_DBD <- matrix(
  0,
  nrow = n_layers,
  ncol = 2 * n_layers + 1
)

for (i in seq_len(n_layers)) {
  
  # Present-day layer DBD
  
  J_layer_DBD[i, i] <-
    (1 - phi_0) / (1 - phi_layer[i])
  
  # Present-day layer porosity
  
  J_layer_DBD[i, n_layers + i] <-
    DBD_layer_depositional[i] /
    (1 - phi_layer[i])
  
  # Depositional reference porosity
  
  J_layer_DBD[i, 2 * n_layers + 1] <-
    -DBD_layer_depositional[i] /
    (1 - phi_0)
}


# ------------------------------------------------------------
# 9. DEPOSITIONAL LAYER DBD COVARIANCE
# ------------------------------------------------------------

Sigma_DBD_layer_depositional <-
  J_layer_DBD %*%
  Sigma_layer_input %*%
  t(J_layer_DBD)

Sigma_DBD_layer_depositional <- (
  Sigma_DBD_layer_depositional +
    t(Sigma_DBD_layer_depositional)
) / 2

uDBD_layer_depositional <- sqrt(
  pmax(diag(Sigma_DBD_layer_depositional), 0)
)


# ------------------------------------------------------------
# 10. CALCULATE CORRECTED SAR
# ------------------------------------------------------------

SAR_corrected <- MAR_original /
  DBD_layer_depositional


# ------------------------------------------------------------
# 11. PROPAGATE SAR UNCERTAINTY
# ------------------------------------------------------------
#
# SAR = MAR / DBD_depositional
#
# Assumes MAR uncertainty is independent of the
# reconstructed depositional DBD uncertainty.
#
# ------------------------------------------------------------

Sigma_MAR <- diag(
  uMAR_original^2,
  nrow = n_layers,
  ncol = n_layers
)

Sigma_SAR_input <- matrix(
  0,
  nrow = 2 * n_layers,
  ncol = 2 * n_layers
)

Sigma_SAR_input[
  1:n_layers,
  1:n_layers
] <- Sigma_MAR

Sigma_SAR_input[
  (n_layers + 1):(2 * n_layers),
  (n_layers + 1):(2 * n_layers)
] <- Sigma_DBD_layer_depositional


# SAR Jacobian

J_SAR <- matrix(
  0,
  nrow = n_layers,
  ncol = 2 * n_layers
)

for (i in seq_len(n_layers)) {
  
  # Derivative with respect to MAR
  
  J_SAR[i, i] <-
    1 / DBD_layer_depositional[i]
  
  # Derivative with respect to depositional DBD
  
  J_SAR[i, n_layers + i] <-
    -MAR_original[i] /
    DBD_layer_depositional[i]^2
}


# Full SAR covariance

Sigma_SAR_corrected <-
  J_SAR %*%
  Sigma_SAR_input %*%
  t(J_SAR)

Sigma_SAR_corrected <- (
  Sigma_SAR_corrected +
    t(Sigma_SAR_corrected)
) / 2

uSAR_corrected <- sqrt(
  pmax(diag(Sigma_SAR_corrected), 0)
)


# ------------------------------------------------------------
# 12. GENERATE INDEPENDENT CF/CRS CORRECTION TABLE
# ------------------------------------------------------------

CF_compaction_table <- data.frame(
  
  `Top of Interval z(i) cm` =
    layer_depths,
  
  `Mid Depth zi cm` =
    as.numeric(CF_table[layer_rows, 2]),
  
  `Age t(i) yr` =
    Age_original,
  
  `Year CE` =
    Year_original,
  
  `Mass Accumulation Rate MAR g/(cm^2 yr)` =
    MAR_original,
  
  `MAR Uncertainty g/(cm^2 yr)` =
    uMAR_original,
  
  `Present Day Layer DBD g/cm^3` =
    DBD_layer_original,
  
  `Present Day Layer DBD Uncertainty g/cm^3` =
    uDBD_layer_present,
  
  `Present Day Layer Porosity` =
    phi_layer,
  
  `Present Day Layer Porosity Uncertainty` =
    u_phi_layer,
  
  `Depositional Layer DBD g/cm^3` =
    DBD_layer_depositional,
  
  `Depositional Layer DBD Uncertainty g/cm^3` =
    uDBD_layer_depositional,
  
  `Original CF SAR cm/yr` =
    SAR_original,
  
  `Original CF SAR Uncertainty cm/yr` =
    uSAR_original,
  
  `Corrected CF SAR cm/yr` =
    SAR_corrected,
  
  `Corrected CF SAR Uncertainty cm/yr` =
    uSAR_corrected,
  
  `SAR Correction Factor` =
    DBD_layer_original / DBD_layer_depositional,
  
  check.names = FALSE
)


# ------------------------------------------------------------
# 13. STORE COVARIANCE MATRICES AS ATTRIBUTES
# ------------------------------------------------------------

attr(
  CF_compaction_table,
  "Layer_Porosity_Interpolation_Weights"
) <- W_phi_layer

attr(
  CF_compaction_table,
  "Layer_Porosity_Covariance"
) <- Sigma_phi_layer

attr(
  CF_compaction_table,
  "Layer_DBD_Interpolation_Weights"
) <- W_DBD_layer

attr(
  CF_compaction_table,
  "Present_Layer_DBD_Covariance"
) <- Sigma_DBD_layer_present

attr(
  CF_compaction_table,
  "Layer_Input_Covariance"
) <- Sigma_layer_input

attr(
  CF_compaction_table,
  "Depositional_Layer_DBD_Jacobian"
) <- J_layer_DBD

attr(
  CF_compaction_table,
  "Depositional_Layer_DBD_Covariance"
) <- Sigma_DBD_layer_depositional

attr(
  CF_compaction_table,
  "Corrected_SAR_Jacobian"
) <- J_SAR

attr(
  CF_compaction_table,
  "Corrected_SAR_Covariance"
) <- Sigma_SAR_corrected


# ------------------------------------------------------------
# 14. EXPORT RESULTS
# ------------------------------------------------------------

folder_name_compaction <- "CF CRS Compaction Correction"

dir.create(
  folder_name_compaction,
  showWarnings = FALSE,
  recursive = TRUE
)

write_xlsx(
  CF_compaction_table,
  file.path(
    folder_name_compaction,
    "CF_CRS_compaction_corrected.xlsx"
  )
)


# Export covariance matrices

CF_covariance_output <- list(
  
  `Layer Porosity` =
    as.data.frame(Sigma_phi_layer),
  
  `Present Layer DBD` =
    as.data.frame(Sigma_DBD_layer_present),
  
  `Depositional Layer DBD` =
    as.data.frame(Sigma_DBD_layer_depositional),
  
  `Corrected SAR` =
    as.data.frame(Sigma_SAR_corrected),
  
  `Porosity Interpolation Weights` =
    as.data.frame(W_phi_layer),
  
  `DBD Interpolation Weights` =
    as.data.frame(W_DBD_layer)
)

write_xlsx(
  CF_covariance_output,
  file.path(
    folder_name_compaction,
    "CF_CRS_compaction_covariance.xlsx"
  )
)
# Depositional SAR plot

pdf(
  file.path(
    folder_name_compaction,
    "Depositional_SAR_plot.pdf"
  ),
  width = 8,
  height = 6
)

SAR_depositional_plot_table <- data.frame(x = CF_compaction_table[["Corrected CF SAR cm/yr"]],x_uncertainty = CF_compaction_table[["Corrected CF SAR Uncertainty cm/yr"]],y = CF_compaction_table[["Year CE"]])

SAR_depositional_plot_table <- na.omit(
  SAR_depositional_plot_table
)

SAR_depositional_plot_table <- SAR_depositional_plot_table[
  order(SAR_depositional_plot_table$y),
]

SAR_depositional_plot <- ggplot() +
  
  geom_point(
    data = SAR_depositional_plot_table,
    aes(x = x, y = y)
  ) +
  
  geom_path(
    data = SAR_depositional_plot_table,
    aes(x = x, y = y),
    color = "blue"
  ) +
  
  geom_errorbar(
    data = SAR_depositional_plot_table,
    aes(
      y = y,
      xmin = x - x_uncertainty,
      xmax = x + x_uncertainty
    ),
    width = 0.2
  ) +
  
  labs(
    title = "Depositional Sediment Accumulation Rate by Year",
    x = "Depositional Sediment Accumulation Rate (cm/yr)",
    y = "Year (CE)"
  ) +
  
  dynamic_scale(
    SAR_depositional_plot_table,
    "x"
  ) +
  
  dynamic_scale(
    SAR_depositional_plot_table,
    "y",
    axis = "y"
  ) +
  
  theme_minimal()

print(SAR_depositional_plot)

dev.off()
}

#save the data
folder_name3 = "CF modeling data"
dir.create(folder_name3)
CF_table_path <- file.path(folder_name3, "CF_table.xlsx")
write_xlsx(as.data.frame(CF_table), CF_table_path)
# year plot vs layer depth-z(i)
pdf(file.path(folder_name3, "age_vs_layer_depth.pdf"), width = 8, height = 6)
CF_age_table = array(data = NA, dim = c(length(na.omit(CF_table[, 9])), 3), dimnames = list(NULL, c("x", "x_uncertainty", "y")))
CF_age_table[, 3] = na.omit(CF_table[1:(nrow(CF_table) - 2), 1])
CF_age_table[, 2] = na.omit(CF_table[, 8])
CF_age_table[, 1] = na.omit(CF_table[, 9])
CF_age_table = as.data.frame(CF_age_table)

CF_age_table_plot <- ggplot(
  CF_age_table,
  aes(x = x, y = y)
) +
  geom_point() +
  geom_errorbar(
    aes(
      xmin = x - x_uncertainty,
      xmax = x + x_uncertainty
    ),
    width = 0.2
  ) +
  labs(
    title = "Age vs Layer Depth",
    x = "Year (CE)",
    y = "Layer Depth (cm)"
  ) +
  dynamic_scale(
    CF_age_table,
    "x",
    reverse = TRUE
  ) +
  dynamic_scale(
    CF_age_table,
    "y",
    axis = "y",
    reverse = TRUE
  ) +
  theme_minimal()

print(CF_age_table_plot)
dev.off()
#age vs mass depth plot
# year plot vs layer depth-z(i)
pdf(file.path(folder_name3, "CF_age_vs_mass_depth.pdf"), width = 8, height = 6)
CF_age_md_table = array(data = NA, dim = c(length(na.omit(CF_table[, 9])), 3), dimnames = list(NULL, c("x", "x_uncertainty", "y")))
CF_age_md_table[, 3] = na.omit(mass_table[1:(nrow(CF_table) - 2), 8])
CF_age_md_table[, 2] = na.omit(CF_table[, 8])
CF_age_md_table[, 1] = na.omit(CF_table[, 9])
CF_age_md_table = as.data.frame(CF_age_md_table)

CF_age_md_table_plot <- ggplot(
  CF_age_md_table,
  aes(x = x, y = y)
) +
  geom_point() +
  geom_errorbar(
    aes(
      xmin = x - x_uncertainty,
      xmax = x + x_uncertainty
    ),
    width = 0.2
  ) +
  labs(
    title = "Age vs Cumulative Mass Depth",
    x = "Year (CE)",
    y = "Cumulative Mass Depth (g/cm^2)"
  ) +
  dynamic_scale(
    CF_age_md_table,
    "x",
    reverse = TRUE
  ) +
  dynamic_scale(
    CF_age_md_table,
    "y",
    axis = "y",
    reverse = TRUE
  ) +
  theme_minimal()

print(CF_age_md_table_plot)
dev.off()

# MAR plot
pdf(file.path(folder_name3, "MAR_plot.pdf"), width = 8, height = 6)
MAR_plot_table = array(data = NA, dim = c(length(na.omit(CF_table[, 9])), 3), dimnames = list(NULL, c("x", "x_uncertainty", "y")))
MAR_plot_table[, 3] = na.omit(CF_table[, 9])
MAR_plot_table[, 2] = na.omit(CF_table[, 11])
MAR_plot_table[, 1] = na.omit(CF_table[, 10])
MAR_plot_table = as.data.frame(MAR_plot_table)

MAR_plot_table <- MAR_plot_table[order(MAR_plot_table$y), ]

MAR_plot <- ggplot() +
  geom_point(
    data = MAR_plot_table,
    aes(x = x, y = y)
  ) +
  geom_path(
    data = MAR_plot_table,
    aes(x = x, y = y),
    color = "red"
  ) +
  geom_errorbar(
    data = MAR_plot_table,
    aes(
      y = y,
      xmin = x - x_uncertainty,
      xmax = x + x_uncertainty
    ),
    width = 0.2
  ) +
  labs(
    title = "Mass Accumulation Rate by Year",
    x = "Mass Accumualtion Rate (g/(cm^2 yr))",
    y = "Year (CE)"
  ) +
  dynamic_scale(
    MAR_plot_table,
    "x"
  ) +
  dynamic_scale(
    MAR_plot_table,
    "y",
    axis = "y"
  ) +
  theme_minimal()

print(MAR_plot)
dev.off()
#SAR plot
pdf(file.path(folder_name3, "SAR_plot.pdf"), width = 8, height = 6)
SAR_plot_table = array(data = NA, dim = c(length(na.omit(CF_table[, 9])), 3), dimnames = list(NULL, c("x", "x_uncertainty", "y")))
SAR_plot_table[, 3] = na.omit(CF_table[, 9])
SAR_plot_table[, 2] = na.omit(CF_table[, 17])
SAR_plot_table[, 1] = na.omit(CF_table[, 16])
SAR_plot_table = as.data.frame(SAR_plot_table)

SAR_plot <- ggplot() +
  geom_point(
    data = SAR_plot_table,
    aes(x = x, y = y)
  ) +
  geom_path(
    data = SAR_plot_table,
    aes(x = x, y = y),
    color = "blue"
  ) +
  geom_errorbar(
    data = SAR_plot_table,
    aes(
      y = y,
      xmin = x - x_uncertainty,
      xmax = x + x_uncertainty
    ),
    width = 0.2
  ) +
  labs(
    title = "Sediment Accumulation Rate by Year",
    x = "Sediment Accumulation Rate (cm/yr)",
    y = "Year (CE)"
  ) +
  dynamic_scale(
    SAR_plot_table,
    "x"
  ) +
  dynamic_scale(
    SAR_plot_table,
    "y",
    axis = "y"
  ) +
  theme_minimal()

print(SAR_plot)
dev.off()
if (exists("model_summaries_surf_DBD")){
write_xlsx(model_summaries_surf_DBD, path = "CF modeling data/DBD_model_summaries.xlsx")
}
# run some CF/CFCS comparisons, save plots and data to their own folder
if (surface_active_zone_check==FALSE){
folder_name5 = "CF_CFCS model comparisons"
dir.create(folder_name5)
CF_CFCS_table=array(data=NA, dim=c(nrow(mass_table),6), dimnames=list(NULL,c("Top of Interval z(i) cm", "Mid Depth zi cm","Mass Depth Section mi gcm^-2","Cumulative Mass Depth Layer m(i) gcm^-2","CF Age CE","CFCS Age CE")))
CF_CFCS_table[,1:2]=mass_table[,1:2]
CF_CFCS_table[,3]=mass_table[,10]
CF_CFCS_table[,4]=mass_table[,8]
CF_CFCS_table[,5]=CF_table[,9]
# get CFCS years into proper format
CFCS_format_fix = CFCS_table[, 5]
if (!is.na(CFCS_format_fix[1])) {
  num_rows_CFCS <- length(CFCS_format_fix)  # Get the number of rows in the vector
  # Initialize a new vector to store the result
  new_CFCS_format_fix <- vector("numeric", length = num_rows_CFCS * 2 + 1)  # Create a numeric vector of double the length plus one for the initial NA
  # Fill the new vector with NA and data values
  new_CFCS_format_fix[1] <- NA  # Add an NA before the first data entry
  for (i in 1:num_rows_CFCS) {
    new_CFCS_format_fix[i * 2] <- CFCS_format_fix[i]  # Fill the data entry
    new_CFCS_format_fix[i * 2 + 1] <- NA  # Add an NA after each data entry
  }
  # Remove the last NA if the original data length is odd
  if (num_rows_CFCS %% 2 == 1) {
    new_CFCS_format_fix <- new_CFCS_format_fix[-length(new_CFCS_format_fix)]
  }
  CFCS_format_fix = new_CFCS_format_fix  # Assign the new vector back
}
CF_CFCS_table[1:length(CFCS_format_fix),6]=CFCS_format_fix
#make the comparison plot
CF_CFCS_table=as.matrix(CF_CFCS_table)
# CF table for depth plot
CF_comp_z_table=na.omit(as.matrix(CF_CFCS_table[, c("Top of Interval z(i) cm", "CF Age CE")]))
#CFCS table for depth plot
CFCS_comp_z_table=na.omit(as.matrix(CF_CFCS_table[, c("Mid Depth zi cm", "CFCS Age CE")]))
#make the plot
pdf(file.path(folder_name5, "age_comp_z.pdf"), width = 8, height = 6)
CF_comp_z_table <- as.data.frame(CF_comp_z_table)
CFCS_comp_z_table <- as.data.frame(CFCS_comp_z_table)

# Rename columns for consistency
colnames(CF_comp_z_table) <- c("Depth", "Age")
colnames(CFCS_comp_z_table) <- c("Depth", "Age")

# Add a new column to identify the source of the data
CF_comp_z_table$Source <- "CF"
CFCS_comp_z_table$Source <- "CFCS"

# Combine the two tables
combined_data <- rbind(CF_comp_z_table, CFCS_comp_z_table)

# Check for missing values
if (any(is.na(combined_data))) {
  print("Missing values found in the data:")
  print(combined_data[is.na(combined_data), ])
}

# Remove rows with missing values
combined_data <- na.omit(combined_data)

# Create the plot with reversed axes
CF_CFCS_comp_plot_z <- ggplot(
  combined_data,
  aes(x = Age, y = Depth, color = Source)
) +
  geom_point() +
  geom_line(
    aes(group = Source)
  ) +
  labs(
    x = "Age (CE)",
    y = "Depth (cm)",
    title = "CF and CFCS Age Model Comparisons (z)"
  ) +
  dynamic_scale(
    combined_data,
    "Age",
    reverse = TRUE
  ) +
  dynamic_scale(
    combined_data,
    "Depth",
    axis = "y",
    reverse = TRUE
  ) +
  scale_color_manual(
    values = c(
      "CF" = "blue",
      "CFCS" = "red"
    )
  ) +
  theme_minimal()

print(CF_CFCS_comp_plot_z)
dev.off()
#mass depth comparison plot
# CF table for mass depth plot
CF_comp_md_table=as.matrix(array(data=NA, dim=c((nrow(concentrations_table)/2),2)))
CF_comp_md_table[,1]=na.omit(mass_table[,10])
CF_comp_md_table[,2]=na.omit(CF_table[,9])
#CFCS table for mass depth plot
CFCS_comp_md_table=as.matrix(array(data=NA, dim=c((nrow(concentrations_table)/2),2)))
CFCS_comp_md_table[,1]=na.omit(mass_table[1:(nrow(mass_table)-2),8])
CFCS_comp_md_table[,2]=na.omit(CFCS_table[,5])
#make the plot
pdf(file.path(folder_name5, "age_comp_md.pdf"), width = 8, height = 6)
CF_comp_md_table <- as.data.frame(CF_comp_md_table)
CFCS_comp_md_table <- as.data.frame(CFCS_comp_md_table)

# Rename columns for consistency
colnames(CF_comp_md_table) <- c("Depth", "Age")
colnames(CFCS_comp_md_table) <- c("Depth", "Age")

# Add a new column to identify the source of the data
CF_comp_md_table$Source <- "CF"
CFCS_comp_md_table$Source <- "CFCS"

# Combine the two tables
combined_data_md <- rbind(CF_comp_md_table, CFCS_comp_md_table)

# Check for missing values
if (any(is.na(combined_data_md))) {
  print("Missing values found in the data:")
  print(combined_data_md[is.na(combined_data_md), ])
}

# Remove rows with missing values
combined_data_md <- na.omit(combined_data_md)

# Create the plot with reversed axes
CF_CFCS_comp_plot_md <- ggplot(
  combined_data_md,
  aes(x = Age, y = Depth, color = Source)
) +
  geom_point() +
  geom_line(
    aes(group = Source)
  ) +
  labs(
    x = "Age (CE)",
    y = "Cumulative Mass Depth (g/cm^2)",
    title = "CF and CFCS Age Model Comparisons (md)"
  ) +
  dynamic_scale(
    combined_data_md,
    "Age",
    reverse = TRUE
  ) +
  dynamic_scale(
    combined_data_md,
    "Depth",
    axis = "y",
    reverse = TRUE
  ) +
  scale_color_manual(
    values = c(
      "CF" = "blue",
      "CFCS" = "red"
    )
  ) +
  theme_minimal()

print(CF_CFCS_comp_plot_md)
dev.off()
# save the cf/cfcs table to CSV
write_xlsx(as.data.frame(CF_CFCS_table), "CF_CFCS model comparisons/CF_CFCS_table.xlsx")
}
#completion dialogue for CF

# Zhang atm flux section
insert_results=c(NA,NA)
est_local_flux <- as.logical(readline(
  prompt = "Do you want to to estimate local atmospheric flux? (TRUE/FALSE): "
))
if (est_local_flux==TRUE){
  flux_table <- data.frame(
    `Lower Latitudinal Bound` = c(70, 60, 50, 40, 30, 20, 10, 0, -20, -30, -40, -50, -70, -80, -90),
    `Upper Latitudinal Bound` = c(80, 70, 60, 50, 40, 30, 20, 10, -10, -20, -30, -40, -60, -70, -80),
    `Equation Slope` = c(0.04483, 0.09515, 0.03139, 0.05889, 0.12995, 0.25347, 0.08746, -0.03217, 0.0603, 0.0133, 0.11289, 0.05475, 0.00955, 0.02246, 0.01689),
    `Equation y Intercept` = c(0.07302, -16.14357, 74.9366, 117.5243, 71.78459, -124.54982, -6.35611, 239.85701, -2.1501, 56.46349, -31.33849, -3.46801, 1.67204, 0.17641, 0.50693),
    `Slope SE`=c(0.00925,0.01825,0.01323,0.01587,0.02924,0.0434,0.01395,0.04825,0.02593,0.02871,0.03322,0.01296,0.00281,0.00543,0.00263),
    `Intercept SE`=c(11.13663,10.59527,17.07211,19.27045,34.97366,97.24497,22.18361,112.95092,34.32702,29.41251,38.30933,17.72279,0.93411,0.97109,0.16005),
    `Pearson r square` = c(0.84, 0.76, 0.25, 0.25, 0.39, 0.59, 0.94, 0.29, 0.66, 0.15, 0.65, 0.80, 0.77, 0.86, 0.92)
  )
  lat_band=as.numeric(readline(prompt = "Enter core latitude (indicate southern latitude with -) "))
  precip_mm=as.numeric(readline(prompt = "Enter local annual precipitation (mm) "))
  for (i in 1:15){
    if (lat_band>= flux_table[i,1] & lat_band<=flux_table[i,2]){
      flux_est=flux_table[i,3]*precip_mm+flux_table[i,4]
      flux_est=(flux_est*60)/10000
      flux_est_unc=sqrt((precip_mm*flux_table[i,5])^2+flux_table[i,6]^2)
      flux_est_unc=(flux_est_unc*60)/10000
      cat("Local curve estimate atmospheric Pb-210 flux is ", flux_est, " +/- ", flux_est_unc, " dpm/(cm^2 y)\n")
      insert_results=c(flux_est,flux_est_unc)
      break
    }
  }
  if ((lat_band < 0 & lat_band > -10) |
      (lat_band < -50 & lat_band > -60)) {
    print("no local curve for atmospheric Pb-210 flux estimate available")
    insert_results=c(NA,NA)
  }
}
finish_CF=FALSE
print("core-derived atmospheric Pb-210 flux (dpm cm^2/yr)")
atm_flux=as.numeric(CF_table[1,5]*decay_const)
atm_flux_uncer=atm_flux*sqrt((as.numeric(CF_table[1,6])/as.numeric(CF_table[1,5]))^2+(decay_const_uncer/decay_const)^2)
print(atm_flux)
flux_table=array(data=NA, dim=c(1,4), dimnames=list(NULL, c("Atmospheric Pb-210 Flux (dpm/(cm^2 yr))", "Atmospheric Pb-210 Flux Uncertainty (dpm/(cm^2 yr))","Estimated Atmospheric Pb-210 Flux (dpm/(cm^2 yr))", "Estimate Pearson's r square")))
flux_table[1,1]=atm_flux
flux_table[1,2]=atm_flux_uncer
flux_table[1,3]=insert_results[1]
flux_table[1,4]=insert_results[2]
flux_table_path <- file.path(folder_name3, "atm_pb210_flux.xlsx")
write_xlsx(as.data.frame(flux_table), flux_table_path)

  print("missing inventory %")
  print(as.numeric((missing_inventory/CF_table[1,5])*100))
  print("surface active zone depth %")
  if (SAZ_depth!=0){
    print(as.numeric((temp_surface_table[SAZ_depth,2]/concentrations_table[nrow(concentrations_table),1])*100))
  }else{
    print(0)
  }
  
  if (as.numeric((missing_inventory/CF_table[1,5])*100)>=20){
    print("missing inventory is greater than 20%! Interpret with caution")
  }
  if (SAZ_depth!=0){
    if (as.numeric((temp_surface_table[SAZ_depth,2]/concentrations_table[nrow(concentrations_table),1])*100)>=20){
      print("surface active zone is greater than 20% of the core! Interpret with caution")
  }
  }
  if (inventory_completed==FALSE){
    print("Could not complete inventory! CF is not a valid model")
  }
  finish_CF=readline(prompt = "CF modeling complete. Continue to CA? (True or False) ")


# save the core boundary issues (SAZ and missing inventory)
core_boundary_info=array(data=NA, dim=c(1,4), dimnames=list(NULL, c("SAZ depth (cm)", "SAZ core %", "missing inventory (dpm/cm^2)", "missing inventory %")))
if (SAZ_depth!=0){
  core_boundary_info[1,1]=temp_surface_table[SAZ_depth,2]
  core_boundary_info[1,2]=as.numeric((temp_surface_table[SAZ_depth,2]/concentrations_table[nrow(concentrations_table),1])*100)
}else{
  core_boundary_info[1,1]=0
  core_boundary_info[1,2]=0
}
core_boundary_info[1,3]=missing_inventory
core_boundary_info[1,4]=(missing_inventory/CF_table[1,5])*100
boundary_table_path <- file.path(folder_name3, "core_boundary_info.xlsx")
write_xlsx(as.data.frame(core_boundary_info), boundary_table_path)
}else{
  finish_CF=readline(prompt = "CF modeling skipped. Do you want to continue to CA? (True or False) ")
}

#start of CA section
if (finish_CF==TRUE){
#initialize array
folder_name4 = "CA modeling data"
dir.create(folder_name4)
CA_table=array(data=NA, dim=c((nrow(concentrations_table)+1),20), dimnames=list(NULL,c("Mid Depth zi cm", "ln(Ci)","Age ti yr","Age Uncertainty u(ti) yr","CA Age CE","ur(l)","ur(C0)","ur(Ci)","t(i) yr","u(t(i)) yr","delta t yr","Uncertainty delta t yr"," Sediment Accumulation Rate s cm/yr","Sediment Accumulation Rate Uncertainty u(s) cm/yr","Mass Accumualtion Rate r g/(cm^2 yr)","Mass Accumualtion Rate Uncertainty u(r) g/(cm^2 yr)","Mass Depth Section mi gcm^-2","Mass Depth Section Uncertainty u(mi) gcm^-2","Cumulative Mass Depth Layer m(i) gcm^-2","Cumulative Mass Depth Layer Uncertainty u(m(i)) gcm^-2")))
colnames(CA_table) <- gsub("delta", "\u0394", colnames(CA_table))
# load copied data
CA_table[,1]=mass_table[,2]
CA_table[1:nrow(concentrations_table),2]=concentrations_table[,13]
CA_table[,17]=mass_table[,10]
CA_table[,18]=mass_table[,11]
CA_table[,19]=mass_table[,8]
CA_table[,20]=mass_table[,9]


#start the loops
#age ti yr, CA year, ur(l), ur(C0), ur(Ci)
for (i in 1:nrow(CA_table)){
  if (i%%2==0){
    CA_table[i,3]=log(concentrations_table[1,8]/concentrations_table[i,6])/decay_const
    CA_table[i,5]=numeric_date-CA_table[i,3]
    CA_table[i,6]=decay_const_uncer*CA_table[i,3]
    CA_table[i,7]=concentrations_table[1,9]/concentrations_table[1,8]
    CA_table[i,8]=concentrations_table[i,7]/concentrations_table[i,6]
  }
}
# age ti uncertainty
for (i in 1:nrow(CA_table)){
  if (i%%2==0){
    CA_table[i,4]=sqrt((CA_table[i,6])^2+(CA_table[i,7])^2+(CA_table[i,8])^2)/decay_const
  }
}
#t(i) yr + uncertainty loop
for (i in 1:(nrow(CA_table))){
  if (i%%2==1){
    if (i==1){
      CA_table[i,9]=0
      CA_table[i,10]=0
    }else if (i==nrow(CA_table)){
      fraction=(CA_table[i,"Cumulative Mass Depth Layer m(i) gcm^-2"]-CA_table[(i-1),"Mass Depth Section mi gcm^-2"])/(CA_table[(i-1),"Mass Depth Section mi gcm^-2"]-CA_table[(i-3),"Mass Depth Section mi gcm^-2"])
      CA_table[i,9]=CA_table[(i-1),"Age ti yr"]+(CA_table[(i-1),"Age ti yr"]-CA_table[(i-3),"Age ti yr"])*fraction
      CA_table[i,10]=sqrt((1+fraction)^2*CA_table[(i-1),"Age Uncertainty u(ti) yr"]^2+fraction^2*CA_table[(i-3),"Age Uncertainty u(ti) yr"]^2)
    }else{
      fraction =
        (CA_table[i,"Cumulative Mass Depth Layer m(i) gcm^-2"] -
           CA_table[(i-1),"Mass Depth Section mi gcm^-2"]) /
        (CA_table[(i+1),"Mass Depth Section mi gcm^-2"] -
           CA_table[(i-1),"Mass Depth Section mi gcm^-2"])
      
      CA_table[i,"t(i) yr"] =
        CA_table[(i-1),"Age ti yr"] +
        ((CA_table[(i+1),"Age ti yr"] -
            CA_table[(i-1),"Age ti yr"])) * fraction
      
      CA_table[i,"u(t(i)) yr"] =
        sqrt((1-fraction)^2 * CA_table[(i-1),"Age Uncertainty u(ti) yr"]^2 +
               (fraction)^2 * CA_table[(i+1),"Age Uncertainty u(ti) yr"]^2)
    }
  }
}
# next set variables loop, delta t yr, uncertainty delta t yr
for (i in 1:(nrow(CA_table))){
  if (i%%2==0){
    CA_table[i,11]=CA_table[(i+1),9]-CA_table[(i-1),9]
    CA_table[i,12]= sqrt((CA_table[(i+1),10])^2+(CA_table[(i-1),10])^2)
  }
}
#SAR, SAR uncertainty, MAR, MAR uncertainty loop
for (i in 1:(nrow(CA_table))){
  if (i%%2==0){
    CA_table[i,13]=(mass_table[(i+1),1]-mass_table[(i-1),1])/CA_table[i,11]
    CA_table[i,14]=(CA_table[i,13]*CA_table[i,12])/CA_table[i,11]
    CA_table[i,15]=(mass_table[(i+1),8]-mass_table[(i-1),8])/CA_table[i,11]
    CA_table[i,16]=(CA_table[i,15]*CA_table[i,12])/CA_table[i,11]
  }
}
if (exists("zi_correction_table")) {
  
  # Output table
  ccSAR_CA <- data.frame(
    `CA Uncompacted SAR (cm/yr)` =
      rep(NA_real_, nrow(CA_table)),
    `CA Uncompacted SAR Uncertainty (cm/yr)` =
      rep(NA_real_, nrow(CA_table))
  )
  
  # Covariance matrix of corrected top-of-interval depths
  Sigma_top <- attr(
    zi_correction_table,
    "Corrected_TopDepth_Covariance"
  )
  
  # Covariance matrix of corrected section thicknesses
  Sigma_hprime <- attr(
    zi_correction_table,
    "Corrected_Thickness_Covariance"
  )
  
  if (is.null(Sigma_top) || is.null(Sigma_hprime)) {
    stop("Corrected-depth covariance matrices are missing.")
  }
  
  n_sections <- nrow(Sigma_hprime)
  
  if (nrow(CA_table) < 2 * n_sections) {
    stop("CA_table has fewer rows than required for the corrected sections.")
  }
  
  # Calculate SAR for each even-numbered CA row
  for (i in seq_len(nrow(CA_table))) {
    
    if (i %% 2 != 0) next
    
    # Corresponding sediment section
    section <- i / 2
    
    if (section > n_sections) next
    
    # Corrected thickness of this section
    delta_z <- as.numeric(
      zi_correction_table[
        section * 2,
        "Uncompacted Section Thickness (cm)"
      ]
    )
    
    # CA age interval and its uncertainty
    delta_t <- as.numeric(CA_table[i, 11])
    sigma_delta_t <- as.numeric(CA_table[i, 12])
    
    # Validate inputs
    if (!is.finite(delta_z) ||
        !is.finite(delta_t) ||
        delta_t == 0 ||
        !is.finite(sigma_delta_t) ||
        sigma_delta_t < 0) {
      next
    }
    
    # Compaction-corrected SAR
    ccSAR_CA[i, 1] <- delta_z / delta_t
    
    # Variance of corrected section thickness
    var_delta_z <- Sigma_hprime[section, section]
    
    # First-order uncertainty propagation:
    # SAR = delta_z / delta_t
    var_SAR <-
      var_delta_z / delta_t^2 +
      (delta_z^2 * sigma_delta_t^2) / delta_t^4
    
    ccSAR_CA[i, 2] <- sqrt(max(var_SAR, 0))
  }
  
  # Export only rows corresponding to CA sediment intervals
  valid_rows <- which(
    seq_len(nrow(CA_table)) %% 2 == 0 &
      is.finite(ccSAR_CA[, 1])
  )
  zi_exp=na.omit(zi_correction_table[,10])
  ccSAR_CA_export <- data.frame(
    zi_exp,
    ccSAR_CA[valid_rows, 1],
    ccSAR_CA[valid_rows, 2],
    check.names = FALSE
  )
  
  names(ccSAR_CA_export) <- c(
    "Uncompacted Mid-Depth (cm)",
    "CA Uncompacted SAR (cm/yr)",
    "CA Uncompacted SAR Uncertainty (cm/yr)"
  )
  
  ccSAR_CA_table_path <- file.path(
    folder_name4,
    "CA_ccSAR.xlsx"
  )
  
  write_xlsx(ccSAR_CA_export, ccSAR_CA_table_path)
  # Plot CA compaction-corrected SAR
  pdf(
    file.path(folder_name4, "ccSAR_plot_CA.pdf"),
    width = 8,
    height = 6
  )
  
  valid_plot_rows <- which(
    is.finite(CA_table[, 5]) &
      is.finite(ccSAR_CA[, 1]) &
      is.finite(ccSAR_CA[, 2])
  )
  
  ccSAR_plot_table_CA <- data.frame(
    x = as.numeric(CA_table[valid_plot_rows, 5]),
    y = as.numeric(ccSAR_CA[valid_plot_rows, 1]),
    y_uncertainty = as.numeric(ccSAR_CA[valid_plot_rows, 2])
  )
  
  ccSAR_plot_CA <- ggplot(
    ccSAR_plot_table_CA,
    aes(x = x, y = y)
  ) +
    geom_point() +
    geom_path(color = "blue") +
    geom_errorbar(
      aes(
        ymin = y - y_uncertainty,
        ymax = y + y_uncertainty
      ),
      width = 0.2
    ) +
    labs(
      title = "Compaction Corrected Sediment Accumulation Rate by Year",
      x = "Year (CE)",
      y = "Sediment Accumulation Rate (cm/yr)"
    ) +
    dynamic_scale(
      ccSAR_plot_table_CA,
      "x",
      reverse = TRUE
    ) +
    dynamic_scale(
      ccSAR_plot_table_CA,
      "y",
      axis = "y"
    ) +
    theme_minimal()
  
  print(ccSAR_plot_CA)
  dev.off()
}
# save off the data
CA_table_path <- file.path(folder_name4, "CA_table.xlsx")
write_xlsx(as.data.frame(CA_table), CA_table_path)

#make the plots
#log actvity/mass depth
pdf(file.path(folder_name4, "excess_log_pb210_md_CA.pdf"), width = 8, height = 6)
elogpb210_md_plotting_table_CA=array(data=NA, dim=c(nrow(CA_table)/2,3), dimnames=list(NULL, c("x","y", "y_uncertainty")))
elogpb210_md_plotting_table_CA[,1]=na.omit(mass_table[,10])
elogpb210_md_plotting_table_CA[,2]=concentrations_table[seq(2, nrow(concentrations_table), by=2),13]
elogpb210_md_plotting_table_CA[,3]=concentrations_table[seq(2, nrow(CA_table), by=2),14]
elogpb210_md_plotting_table_CA=as.data.frame(elogpb210_md_plotting_table_CA)
elogpb210_md_plot_CA <- ggplot(
  elogpb210_md_plotting_table_CA,
  aes(x = x, y = y)
) +
  geom_point() +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    color = "red"
  ) +
  labs(
    title = "Excess Log-scale Pb-210 vs Mass Depth",
    x = "Cumulative Mass Depth (g/cm^2)",
    y = "Excess Log-scale Pb-210 (dpm/g)"
  ) +
  dynamic_scale(
    elogpb210_md_plotting_table_CA,
    "x"
  ) +
  dynamic_scale(
    elogpb210_md_plotting_table_CA,
    "y",
    axis = "y"
  ) +
  theme_minimal()

print(elogpb210_md_plot_CA)
dev.off()
#log actvity/depth
pdf(file.path(folder_name4, "excess_log_pb210_CA.pdf"), width = 8, height = 6)
elogpb210_plotting_table_CA=array(data=NA, dim=c(nrow(CA_table)/2,3), dimnames=list(NULL, c("x","y", "y_uncertainty")))
elogpb210_plotting_table_CA[,1]=na.omit(CA_table[,1])
elogpb210_plotting_table_CA[,2]=concentrations_table[seq(2, nrow(CA_table), by=2),13]
elogpb210_plotting_table_CA[,3]=concentrations_table[seq(2, nrow(CA_table), by=2),14]
elogpb210_plotting_table_CA=as.data.frame(elogpb210_plotting_table_CA)
elogpb210_CA_plot <- ggplot(
  elogpb210_plotting_table_CA,
  aes(x = x, y = y)
) +
  geom_point() +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    color = "red"
  ) +
  labs(
    title = "Excess Log-scale Pb-210 vs Depth",
    x = "Mid Depth (cm)",
    y = "Excess Log-scale Pb-210 (dpm/g)"
  ) +
  dynamic_scale(
    elogpb210_plotting_table_CA,
    "x"
  ) +
  dynamic_scale(
    elogpb210_plotting_table_CA,
    "y",
    axis = "y"
  ) +
  theme_minimal()

print(elogpb210_CA_plot)
dev.off()


## make CA SAR plot
pdf(file.path(folder_name4, "SAR_plot_CA.pdf"), width = 8, height = 6)
SAR_plot_table_CA = array(data = NA, dim = c((length(na.omit(CA_table[, 13]))), 3), dimnames = list(NULL, c("x", "y", "y_uncertainty")))
SAR_plot_table_CA[, 1] = na.omit(CA_table[1:(length(CA_table[, 13])), 5])
SAR_plot_table_CA[, 2] = na.omit(CA_table[1:(length(CA_table[, 13])), 13])
SAR_plot_table_CA[, 3] = na.omit(CA_table[1:(length(CA_table[, 13])), 14])
SAR_plot_table_CA = as.data.frame(SAR_plot_table_CA)

SAR_plot_CA <- ggplot() +
  geom_point(
    data = SAR_plot_table_CA,
    aes(x = x, y = y)
  ) +
  geom_path(
    data = SAR_plot_table_CA,
    aes(x = x, y = y),
    color = "blue"
  ) +
  geom_errorbar(
    data = SAR_plot_table_CA,
    aes(
      x = x,
      ymin = y - y_uncertainty,
      ymax = y + y_uncertainty
    ),
    width = 0.2
  ) +
  labs(
    title = "Sediment Accumulation Rate by Year",
    x = "Year (CE)",
    y = "Sediment Accumulation Rate cm/yr"
  ) +
  dynamic_scale(
    SAR_plot_table_CA,
    "x",
    reverse = TRUE
  ) +
  dynamic_scale(
    SAR_plot_table_CA,
    "y",
    axis = "y"
  ) +
  theme_minimal()

print(SAR_plot_CA)
dev.off()

## make CA MAR plot
pdf(file.path(folder_name4, "MAR_plot_CA.pdf"), width = 8, height = 6)
MAR_plot_table_CA = array(data = NA, dim = c((length(na.omit(CA_table[, 13]))), 3), dimnames = list(NULL, c("x", "y", "y_uncertainty")))
MAR_plot_table_CA[, 1] = na.omit(CA_table[1:(length(CA_table[, 13])), 5])
MAR_plot_table_CA[, 2] = na.omit(CA_table[1:(length(CA_table[, 13])), 15])
MAR_plot_table_CA[, 3] = na.omit(CA_table[1:(length(CA_table[, 13])), 16])
MAR_plot_table_CA = as.data.frame(MAR_plot_table_CA)

MAR_plot_CA <- ggplot() +
  geom_point(
    data = MAR_plot_table_CA,
    aes(x = x, y = y)
  ) +
  geom_path(
    data = MAR_plot_table_CA,
    aes(x = x, y = y),
    color = "red"
  ) +
  geom_errorbar(
    data = MAR_plot_table_CA,
    aes(
      x = x,
      ymin = y - y_uncertainty,
      ymax = y + y_uncertainty
    ),
    width = 0.2
  ) +
  labs(
    title = "Mass Accumulation Rate by Year",
    x = "Year (CE)",
    y = "Mass Accumulation Rate (g/(cm^2 yr))"
  ) +
  dynamic_scale(
    MAR_plot_table_CA,
    "x",
    reverse = TRUE
  ) +
  dynamic_scale(
    MAR_plot_table_CA,
    "y",
    axis = "y"
  ) +
  theme_minimal()

print(MAR_plot_CA)
dev.off()
#age-depth plot
pdf(file.path(folder_name4, "age_vs_mid_depth.pdf"), width = 8, height = 6)
CA_age_table = array(data = NA, dim = c(length(na.omit(CA_table[, 5])), 3), dimnames = list(NULL, c("y", "x_uncertainty", "x")))
CA_age_table[, 1] = na.omit(CA_table[, 1])
CA_age_table[, 2] = na.omit(CA_table[, 4])
CA_age_table[, 3] = na.omit(CA_table[, 5])
CA_age_table = as.data.frame(CA_age_table)

CA_age_table_plot <- ggplot(
  CA_age_table,
  aes(x = x, y = y)
) +
  geom_point() +
  geom_errorbar(
    aes(
      xmin = x - x_uncertainty,
      xmax = x + x_uncertainty
    ),
    width = 0.2
  ) +
  labs(
    title = "Age vs Mid Depth",
    x = "Year (CE)",
    y = "Mid Depth (cm)"
  ) +
  dynamic_scale(
    CA_age_table,
    "x",
    reverse = TRUE
  ) +
  dynamic_scale(
    CA_age_table,
    "y",
    axis = "y",
    reverse = TRUE
  ) +
  theme_minimal()

print(CA_age_table_plot)
dev.off()

#age mass-depth plot
pdf(file.path(folder_name4, "age_vs_mass_depth.pdf"), width = 8, height = 6)
CA_age_md_table = array(data = NA, dim = c(length(na.omit(CA_table[, 5])), 3), dimnames = list(NULL, c("y", "x_uncertainty", "x")))
CA_age_md_table[, 1] = na.omit(CA_table[, 17])
CA_age_md_table[, 2] = na.omit(CA_table[, 4])
CA_age_md_table[, 3] = na.omit(CA_table[, 5])
CA_age_md_table = as.data.frame(CA_age_md_table)

CA_age_md_table_plot <- ggplot(
  CA_age_md_table,
  aes(x = x, y = y)
) +
  geom_point() +
  geom_errorbar(
    aes(
      xmin = x - x_uncertainty,
      xmax = x + x_uncertainty
    ),
    width = 0.2
  ) +
  labs(
    title = "Age vs Mass Depth",
    x = "Year (CE)",
    y = "Cumulative Mass Depth (g/cm^2)"
  ) +
  dynamic_scale(
    CA_age_md_table,
    "x",
    reverse = TRUE
  ) +
  dynamic_scale(
    CA_age_md_table,
    "y",
    axis = "y",
    reverse = TRUE
  ) +
  theme_minimal()

print(CA_age_md_table_plot)
dev.off()

#CA model validity checker
valid_check=0
for (i in 1:nrow(CA_table)){
  if (i%%2==0){
    valid_check=CA_age_table[i,3]-CA_age_table[(i+2),3]
    if(valid_check<0){
      print("Age-inversions detected, CA model assumptions violated! CA is not a valid model")
      break
    }
  }
}
}
#optional accumulation rate calculator for other materials in a core
material_filepaths=c()
material_accum_opt=readline(prompt = "CA modeling complete. Do you wish to calculate accumulation rates for other materials present in the core? (True or False) ")
mat_accum_finished=TRUE
while (mat_accum_finished==TRUE){
  
  if (material_accum_opt==TRUE){
    print("select material data")
    file_path_new_mat=file.choose()
    material_filepaths=c(material_filepaths,file_path_new_mat)
    new_mat_dat <- read_excel(file_path_new_mat, col_names = TRUE, col_types = "numeric", .name_repair = "unique_quiet")
    new_mat_dat=as.matrix(new_mat_dat)
    new_mat_table=array(NA, dim=c((nrow(new_mat_dat)*2),21))
    first_entry=readline(prompt = "Label your material of interest and provide units (g/g) or (p/g)")
    second_entry <- sub(
      " (\\([^)]*\\))$",
      " Uncertainty \\1",
      first_entry
    )
    third_entry=paste0("Layer ", first_entry)
    fourth_entry=paste0("Layer ", second_entry)
    ## 1. Get first word (XX)
    xx <- strsplit(first_entry, "\\s+")[[1]][1]
    ## 2. Determine units (ZZ)
    if (grepl("\\(g/g\\)", first_entry)) {
      zz <- "g/(m^2 yr)"
    } else if (grepl("\\(p/g\\)", first_entry)) {
      zz <- "p/(m^2 yr)"
    } else {
      stop("Unknown units in input string.")
    }
    ## 3. Prefixes AAAA
    aaaa <- c("CF", "CFCS", "CA")
    ## 4. Build labels
    acc_labels <- paste(aaaa, xx, "Accumulation Rate", zz)
    acc_labels_unc=paste(aaaa, xx, "Accumulation Rate Uncertainty", zz)
    new_column_names_mat_data=c("Cumulative Mass Depth Layer m(i) gcm^-2","Mass Depth Section mi gcm^-2", first_entry, second_entry,third_entry,fourth_entry, "CF Mass Accumulation g/(m^2 yr)","CF Mass Accumulation Uncertainty g/(m^2 yr)","CF Layer Age (yr CE)",acc_labels[1],acc_labels_unc[1], "CFCS Mass Accumulation g/(m^2 yr)","CFCS Mass Accumulation Uncertainty g/(m^2 yr)","CFCS Section Age (yr CE)",acc_labels[2],acc_labels_unc[2],"CA Mass Accumulation g/(m^2 yr)","CA Mass Accumulation Uncertainty g/(m^2 yr)","CA Section Age (yr CE)", acc_labels[3],acc_labels_unc[3])
    colnames(new_mat_table)=new_column_names_mat_data
    new_mat_table[,"Cumulative Mass Depth Layer m(i) gcm^-2"]=mass_table[1:nrow(new_mat_table),8]
    new_mat_table[,"Mass Depth Section mi gcm^-2"]=mass_table[1:nrow(new_mat_table),10]
    new_mat_table[,"CF Mass Accumulation g/(m^2 yr)"]=CF_table[1:nrow(new_mat_table),10]*1000
    new_mat_table[,"CF Mass Accumulation Uncertainty g/(m^2 yr)"]=CF_table[1:nrow(new_mat_table),11]*1000
    new_mat_table[,"CF Layer Age (yr CE)"]=CF_table[1:nrow(new_mat_table),9]
    #add reformatted data
    for (i in 1:nrow(new_mat_table)){
      if (i %% 2 == 0){
        new_mat_table[i,first_entry]=new_mat_dat[(i/2),1]
        new_mat_table[i,second_entry]=new_mat_dat[(i/2),2]
      }
    }
    original_section_valid <- !is.na(new_mat_table[[first_entry]])
    #now interpolate any missing values
    for (i in 1:nrow(new_mat_table)){
      if (i %% 2 == 0){
        if (is.na(new_mat_table[i,first_entry])){
          
          #find the previous section with a concentration value
          j = i - 2
          while (j >= 2 && !original_section_valid[j]){
            j = j - 2
          }
          
          #find the next section with a concentration value
          k = i + 2
          while (k <= nrow(new_mat_table) && !original_section_valid[k]){
            k = k + 2
          }
          
          #interpolate only if there is a valid section on both sides
          if (j >= 2 && k <= nrow(new_mat_table)){
            
            fraction =
              (new_mat_table[i,"Mass Depth Section mi gcm^-2"] -
                 new_mat_table[j,"Mass Depth Section mi gcm^-2"]) /
              (new_mat_table[k,"Mass Depth Section mi gcm^-2"] -
                 new_mat_table[j,"Mass Depth Section mi gcm^-2"])
            
            new_mat_table[i,first_entry] =
              new_mat_table[j,first_entry] +
              (new_mat_table[k,first_entry] -
                 new_mat_table[j,first_entry]) * fraction
            
            new_mat_table[i,second_entry] =
              sqrt((1-fraction)^2 * new_mat_table[j,second_entry]^2 +
                     fraction^2 * new_mat_table[k,second_entry]^2)
          }
        }
      }
    }
    
    #now interpolate the layer concentrations from the section concentrations
    for (i in 1:nrow(new_mat_table)){
      if (i == 1){
        fraction=new_mat_table[(i+1),"Mass Depth Section mi gcm^-2"]/(new_mat_table[(i+3),"Mass Depth Section mi gcm^-2"]-new_mat_table[(i+1),"Mass Depth Section mi gcm^-2"])
        new_mat_table[i,third_entry]=new_mat_table[(i+1),first_entry]-((new_mat_table[(i+3),first_entry]-new_mat_table[(i+1),first_entry]))*fraction
        new_mat_table[i, fourth_entry]=sqrt((1-fraction)^2*new_mat_table[(i+1),second_entry]^2+(fraction)^2*new_mat_table[(i+3),second_entry]^2)
        
      } else if (i == nrow(new_mat_table)){
        # last row is empty
        
      } else if (i %% 2 == 1){
        fraction =
          (new_mat_table[i,"Cumulative Mass Depth Layer m(i) gcm^-2"] -
             new_mat_table[(i-1),"Mass Depth Section mi gcm^-2"]) /
          (new_mat_table[(i+1),"Mass Depth Section mi gcm^-2"] -
             new_mat_table[(i-1),"Mass Depth Section mi gcm^-2"])
        
        new_mat_table[i,third_entry] =
          new_mat_table[(i-1),first_entry] +
          ((new_mat_table[(i+1),first_entry] -
              new_mat_table[(i-1),first_entry])) * fraction
        
        new_mat_table[i,fourth_entry] =
          sqrt((1-fraction)^2 * new_mat_table[(i-1),second_entry]^2 +
                 (fraction)^2 * new_mat_table[(i+1),second_entry]^2)
        
      }
    }
    
    new_mat_table[,acc_labels[1]]=new_mat_table[,"CF Mass Accumulation g/(m^2 yr)"]*new_mat_table[,third_entry]
    new_mat_table[,acc_labels_unc[1]]=new_mat_table[,acc_labels[1]]*sqrt((new_mat_table[,"CF Mass Accumulation Uncertainty g/(m^2 yr)"]/new_mat_table[,"CF Mass Accumulation g/(m^2 yr)"])^2+(new_mat_table[,fourth_entry]/new_mat_table[,third_entry])^2)
    
    if (surface_active_zone_check==TRUE){
    }else{
      for (i in seq_len(number_of_models)) {
        
        start_j <- if (i == 1) extents[i] * 2 else (extents[i] + 1) * 2
        end_j   <- extents[i + 1] * 2
        
        for (j in seq(start_j, end_j, by = 2)) {
          
          new_mat_table[j, "CFCS Mass Accumulation g/(m^2 yr)"] <-
            MAR_rates[i, 1] * 1000
          
          new_mat_table[j, "CFCS Mass Accumulation Uncertainty g/(m^2 yr)"] <-
            MAR_rates[i, 2] * 1000
        }
      }
      for (k in 1:nrow(new_mat_table)){
        if (k%%2==0){
          new_mat_table[k,"CFCS Section Age (yr CE)"]=CFCS_table[(k/2),5]
        }
      }
      new_mat_table[,acc_labels[2]]=new_mat_table[,"CFCS Mass Accumulation g/(m^2 yr)"]*new_mat_table[,first_entry]
      new_mat_table[,acc_labels_unc[2]]=new_mat_table[,acc_labels[2]]*sqrt((new_mat_table[,"CFCS Mass Accumulation Uncertainty g/(m^2 yr)"]/new_mat_table[,"CFCS Mass Accumulation g/(m^2 yr)"])^2+(new_mat_table[,second_entry]/new_mat_table[,first_entry])^2)
    }
    
    new_mat_table[,"CA Mass Accumulation g/(m^2 yr)"]=CA_table[1:nrow(new_mat_table),15]*1000
    new_mat_table[,"CA Mass Accumulation Uncertainty g/(m^2 yr)"]=CA_table[1:nrow(new_mat_table),16]*1000
    new_mat_table[,"CA Section Age (yr CE)"]=CA_table[1:nrow(new_mat_table),5]
    new_mat_table[,acc_labels[3]]=new_mat_table[,"CA Mass Accumulation g/(m^2 yr)"]*new_mat_table[,first_entry]
    new_mat_table[,acc_labels_unc[3]]=new_mat_table[,acc_labels[3]]*sqrt((new_mat_table[,"CA Mass Accumulation Uncertainty g/(m^2 yr)"]/new_mat_table[,"CA Mass Accumulation g/(m^2 yr)"])^2+(new_mat_table[,second_entry]/new_mat_table[,first_entry])^2)
    
    CA_valid=TRUE
    if(valid_check<0 && material_accum_opt==TRUE){
      CA_valid=FALSE
    }
    main_title <- paste(xx, "Accumulation Rate")
    new_mat_accum_plot_md <- function() {
      
      has_CFCS <- any(!is.na(new_mat_table[, "CFCS Section Age (yr CE)"]))
      has_CA   <- (CA_valid == TRUE)
      
      # CF
      x1 <- new_mat_table[, acc_labels[1]]
      y1 <- new_mat_table[, "Cumulative Mass Depth Layer m(i) gcm^-2"]
      e1 <- new_mat_table[, acc_labels_unc[1]]
      
      # CFCS
      x2 <- new_mat_table[, acc_labels[2]]
      y2 <- new_mat_table[, "Mass Depth Section mi gcm^-2"]
      e2 <- new_mat_table[, acc_labels_unc[2]]
      
      # CA
      if (has_CA) {
        x3 <- new_mat_table[, acc_labels[3]]
        y3 <- new_mat_table[, "Mass Depth Section mi gcm^-2"]
        e3 <- new_mat_table[, acc_labels_unc[3]]
      }
      
      # All y values used in the plot
      y_all <- c(
        y1,
        if (has_CFCS) y2 else NULL,
        if (has_CA) y3 else NULL
      )
      
      # All x values
      x_all <- c(
        x1,
        if (has_CFCS) x2 else NULL,
        if (has_CA) x3 else NULL
      )
      
      # Include x error bars when determining x limits
      x_all_with_error <- c(
        x1 - e1, x1 + e1,
        if (has_CFCS) c(x2 - e2, x2 + e2) else NULL,
        if (has_CA) c(x3 - e3, x3 + e3) else NULL
      )
      
      # Dynamic y-axis
      y_scale <- dynamic_scale(
        data.frame(y = y_all),
        "y",
        axis = "base",
        reverse = TRUE
      )
      
      # Plot limits
      ylim_raw <- range(y_all, na.rm = TRUE)
      xlim <- range(x_all_with_error, na.rm = TRUE)
      
      ylab <- "Cumulative Mass Depth (g/cm^2)"
      xlab <- paste(xx, "Accumulation Rate", zz)
      
      plot(
        x1, y1,
        type = "b",
        pch = 16,
        col = "black",
        xlim = xlim,
        ylim = rev(ylim_raw),
        xlab = xlab,
        ylab = ylab,
        main = main_title,
        yaxt = "n"
      )
      
      # Dynamic reversed y-axis
      axis(
        side = 2,
        at = y_scale$breaks,
        labels = scales::label_number(
          accuracy = y_scale$accuracy
        )(y_scale$breaks),
        las = 2
      )
      
      # CF error bars
      arrows(
        x0 = x1 - e1,
        y0 = y1,
        x1 = x1 + e1,
        y1 = y1,
        angle = 90,
        code = 3,
        length = 0.05,
        col = "black"
      )
      
      # CFCS
      if (has_CFCS) {
        
        points(
          x2, y2,
          pch = 16,
          col = "red"
        )
        
        lines(
          x2, y2,
          col = "red"
        )
        
        arrows(
          x0 = x2 - e2,
          y0 = y2,
          x1 = x2 + e2,
          y1 = y2,
          angle = 90,
          code = 3,
          length = 0.05,
          col = "red"
        )
      }
      
      # CA
      if (has_CA) {
        
        points(
          x3, y3,
          pch = 16,
          col = "blue"
        )
        
        lines(
          x3, y3,
          col = "blue"
        )
        
        arrows(
          x0 = x3 - e3,
          y0 = y3,
          x1 = x3 + e3,
          y1 = y3,
          angle = 90,
          code = 3,
          length = 0.05,
          col = "blue"
        )
      }
      
      # Legend
      legend_labels <- c(acc_labels[1])
      legend_cols   <- c("black")
      
      if (has_CFCS) {
        legend_labels <- c(legend_labels, acc_labels[2])
        legend_cols   <- c(legend_cols, "red")
      }
      
      if (has_CA) {
        legend_labels <- c(legend_labels, acc_labels[3])
        legend_cols   <- c(legend_cols, "blue")
      }
      
      legend(
        "bottomright",
        legend = legend_labels,
        col = legend_cols,
        lty = 1,
        pch = 16,
        bty = "n"
      )
    }
    new_mat_accum_plot_yr <- function() {
      
      has_CFCS <- any(!is.na(new_mat_table[, "CFCS Section Age (yr CE)"]))
      has_CA   <- (CA_valid == TRUE)
      
      # CF
      x1 <- new_mat_table[, acc_labels[1]]
      y1 <- new_mat_table[, "CF Layer Age (yr CE)"]
      e1 <- new_mat_table[, acc_labels_unc[1]]
      
      # CFCS
      x2 <- new_mat_table[, acc_labels[2]]
      y2 <- new_mat_table[, "CFCS Section Age (yr CE)"]
      e2 <- new_mat_table[, acc_labels_unc[2]]
      
      # CA
      if (has_CA) {
        x3 <- new_mat_table[, acc_labels[3]]
        y3 <- new_mat_table[, "CA Section Age (yr CE)"]
        e3 <- new_mat_table[, acc_labels_unc[3]]
      }
      
      # Combine y values for dynamic y-axis calculation
      y_all <- c(
        y1,
        if (has_CFCS) y2 else NULL,
        if (has_CA) y3 else NULL
      )
      
      # Combine x values if needed
      x_all <- c(
        x1,
        if (has_CFCS) x2 else NULL,
        if (has_CA) x3 else NULL
      )
      
      # Include error bars when determining x-axis limits
      x_all_with_error <- c(
        x1 - e1,
        x1 + e1,
        if (has_CFCS) c(x2 - e2, x2 + e2) else NULL,
        if (has_CA) c(x3 - e3, x3 + e3) else NULL
      )
      
      # Determine plot limits
      ylim_raw <- range(y_all, na.rm = TRUE)
      xlim <- range(x_all_with_error, na.rm = TRUE)
      
      # Dynamic y-axis
      y_scale <- dynamic_scale(
        data.frame(y = y_all),
        "y",
        axis = "base"
      )
      
      ylab <- "Year (CE)"
      xlab <- paste(xx, "Accumulation Rate", zz)
      
      # Main CF plot
      plot(
        x1,
        y1,
        type = "b",
        pch = 16,
        col = "black",
        xlim = xlim,
        ylim = ylim_raw,
        xlab = xlab,
        ylab = ylab,
        main = main_title,
        yaxt = "n"
      )
      
      # Dynamic y-axis
      axis(
        side = 2,
        at = y_scale$breaks,
        labels = scales::label_number(
          accuracy = y_scale$accuracy
        )(y_scale$breaks),
        las = 2
      )
      
      # CF error bars
      arrows(
        x0 = x1 - e1,
        y0 = y1,
        x1 = x1 + e1,
        y1 = y1,
        angle = 90,
        code = 3,
        length = 0.05,
        col = "black"
      )
      
      # CFCS
      if (has_CFCS) {
        
        points(
          x2,
          y2,
          pch = 16,
          col = "red"
        )
        
        lines(
          x2,
          y2,
          col = "red"
        )
        
        # CFCS error bars
        arrows(
          x0 = x2 - e2,
          y0 = y2,
          x1 = x2 + e2,
          y1 = y2,
          angle = 90,
          code = 3,
          length = 0.05,
          col = "red"
        )
      }
      
      # CA
      if (has_CA) {
        
        points(
          x3,
          y3,
          pch = 16,
          col = "blue"
        )
        
        lines(
          x3,
          y3,
          col = "blue"
        )
        
        # CA error bars
        arrows(
          x0 = x3 - e3,
          y0 = y3,
          x1 = x3 + e3,
          y1 = y3,
          angle = 90,
          code = 3,
          length = 0.05,
          col = "blue"
        )
      }
      
      # Legend
      legend_labels <- c(acc_labels[1])
      legend_cols   <- c("black")
      
      if (has_CFCS) {
        legend_labels <- c(
          legend_labels,
          acc_labels[2]
        )
        
        legend_cols <- c(
          legend_cols,
          "red"
        )
      }
      
      if (has_CA) {
        legend_labels <- c(
          legend_labels,
          acc_labels[3]
        )
        
        legend_cols <- c(
          legend_cols,
          "blue"
        )
      }
      
      legend(
        "bottomright",
        legend = legend_labels,
        col = legend_cols,
        lty = 1,
        pch = 16,
        bty = "n"
      )
    }
    # build folder name using the naming convention
    folder_name <- paste(xx, "Accumulation Rate")
    
    # full path in current working directory
    folder_path <- file.path(getwd(), folder_name)
    
    # create folder (won't error if it already exists)
    dir.create(folder_path, showWarnings = FALSE)
    
    # 2. To view on screen
    new_mat_accum_plot_yr()
    
    # 3. To save as PDF in working directory
    xlab <- paste(xx, "Accumulation Rate by Year")
    
    # make a safe filename
    plot_file_name <- paste0(gsub(" ", "_", xlab), ".pdf")
    
    # save the plot
    pdf(file.path(folder_path, plot_file_name))
    new_mat_accum_plot_yr()
    dev.off()
    
    
    # 2. To view on screen
    new_mat_accum_plot_md()
    
    # 3. To save as PDF in working directory
    xlab <- paste(xx, "Accumulation Rate by Mass Depth")
    
    # make a safe filename
    plot_file_name <- paste0(gsub(" ", "_", xlab), ".pdf")
    
    # save the plot
    pdf(file.path(folder_path, plot_file_name))
    new_mat_accum_plot_md()
    dev.off()
    
    xlab <- paste(xx, "Accumulation Rate")
    
    # turn that into a safe filename (no spaces)
    file_name <- paste0(gsub(" ", "_", xlab), ".xlsx")
    
    # write new_mat_table to Excel in the working directory
    write_xlsx(
      list(Sheet1 = as.data.frame(new_mat_table)),
      path = file.path(folder_path, file_name)
    )
  }
  mat_accum_finished=readline(prompt = "Accumulation rate calculations complete. Do you wish to calculate accumulation rates for another material (True or False) ")
}


print("modeling complete!")

save.image(file = "pb210analysis_environment.RData")
rm(list = ls(all.names = TRUE))