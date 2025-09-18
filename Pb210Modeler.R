#For updates and bug fixes, please see the actively maintained version of Pb210Modeler:
#https://github.com/ram-24-mp/Pb210Modeler

#packages needed
library(readxl)
library(writexl)
library(ggplot2)
library(scales)
library(segmented)
library(dplyr)
library(zoo)
library(changepoint)

# some useful functions for handling plots dynamically
accuracy_finder <- function (column){
  number_to_check=min_by_lowest_precision(column)
  true_accuracy=round_to_first_sig_digit(number_to_check)
  return (true_accuracy)
}
round_to_first_sig_digit <- function(x) {
  if (x == 0) return(0)
  exponent <- floor(log10(abs(x)))
  lower <- 10^exponent
  upper <- 10^(exponent + 1)
  cutoff <- 5 * lower
  if (x < cutoff) {
    return(lower)
  } else {
    return(upper)
  }
}

min_by_lowest_precision <- function(vec) {
  # Remove NA
  vec <- vec[!is.na(vec)]
  #Function to count digits after decimal
  count_precision <- function(x) {
    if (x %% 1 == 0) return(0)
    nchar(sub("^\\d+\\.", "", sub("0+$", "", as.character(x))))
  }
  #Get minimum decimal precision
  precisions <- sapply(vec, count_precision)
  min_precision <- min(precisions)
  #Define a comparison function based on that precision
  compare_val <- function(x) round(x, digits = min_precision)
  #Find the value in the original vector with the lowest rounded value
  rounded_vec <- sapply(vec, compare_val)
  min_index <- which.min(rounded_vec)
  #Return the original value at that index
  vec[min_index]
}

leading_decimal_zeros <- function(x) {
  if (x == 0) return(0)
  # Ensure fixed-point formatting (avoid scientific notation)
  x_str <- format(abs(x), scientific = FALSE, trim = TRUE)
  #Get everything after the decimal point
  decimal_part <- sub("^\\d+\\.", "", x_str)
  #Find position of first non-zero digit
  match_pos <- regexpr("[1-9]", decimal_part)
  return(as.integer(match_pos))  # strip attributes
}
breaks_by_value <- function(table, column_name) {
  col_data <- table[[column_name]]
  range <- max(col_data, na.rm = TRUE) - min(col_data, na.rm = TRUE)
  step_size <- range / 10
  rounded <- round(step_size, digits = leading_decimal_zeros(step_size))
  return(rounded)
}

# make a numeric date
day=as.numeric(readline(prompt = "Enter Day "))
month=as.numeric(readline(prompt = "Enter Month "))
year=as.numeric(readline(prompt = "Enter Year "))
numeric_date=year+(month-1)/12+(day-1)/365.24

# core dimensions, in cm
dim_defaults=readline(prompt = "Accept Core Dimension Defaults (True or False) ")
if (dim_defaults==TRUE){
  internal_diameter=10.2
  internal_diameter_uncer=0.1
  surface_area=pi*( internal_diameter/2)^2
  surface_area_uncer=pi*internal_diameter*internal_diameter_uncer/2
  interval_thickness=1
} else {
  internal_diameter=as.numeric(readline(prompt = "Enter Diameter (default 10.2) "))
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

# masses
# dry bulk density g cm^-3 & corresponding uncertainty columns are provided by user. Data must be entered with empty row between each entry
mass_table=array(data=NA, dim=c((((core_length*2)/interval_thickness)+1), 11), dimnames=list(NULL, c("Top of Interval z(i) (cm)", "Mid Depth zi (cm)","Interval Thickness delta zi (cm)","Dry Bulk Density DBD gcm^-3","Dry Bulk Density Uncertainty u(DBD) gcm^-3","Aerial Dry Mass delta mi/S gcm^-2","Aerial Dry Mass Uncertainty u(delta mi/S) gcm^-2","Cumulative Mass Depth Layer m(i) gcm^-2","Cumulative Mass Depth Layer Uncertainty u(m(i)) gcm^-2","Mass Depth Section mi gcm^-2","Mass Depth Section Uncertainty u(mi) gcm^-2")))

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

#load in mass information
# Define the path to your Excel file
print("select mass data")
file_path=file.choose()
# Specify the new column names
new_column_names_mass_data <- c("X", "X.1")
# Read the Excel file
mass_data <- read_excel(file_path, col_names = TRUE, col_types = "numeric")
mass_data=as.data.frame(mass_data)
colnames(mass_data)=new_column_names_mass_data

# Perform linear interpolation on the X values
mass_data$interpolated_X <- na.approx(mass_data$X, na.rm = FALSE)

# Initialize a vector for propagated uncertainties
mass_data$propagated_uncertainties <- NA

# Loop through the interpolated X to calculate uncertainties
for (i in 1:nrow(mass_data)) {
  # Only calculate propagated uncertainty if the current uncertainty is NA
  if (is.na(mass_data$X.1[i])) {
    # Find the nearest known X on both sides
    left_index <- i - 1
    while (left_index > 0 && (is.na(mass_data$X[left_index]) || is.na(mass_data$X.1[left_index]))) {
      left_index <- left_index - 1
    }
    
    right_index <- i + 1
    while (right_index < nrow(mass_data) && (is.na(mass_data$X[right_index]) || is.na(mass_data$X.1[right_index]))) {
      right_index <- right_index + 1
    }
    
    # Check if we found valid indices
    if (left_index > 0 && right_index < nrow(mass_data)) {
      delta_y1 <- mass_data$X.1[left_index]  # Use the correct column for uncertainties
      x1 <- left_index
      delta_y2 <- mass_data$X.1[right_index]  # Use the correct column for uncertainties
      x2 <- right_index
      partial_der=(i-x1)/(x2-x1)
      
      # Calculate propagated uncertainty
      propagated_uncertainty <- sqrt(((1-partial_der)*delta_y1)^2 + (partial_der * delta_y2)^2)
      mass_data$propagated_uncertainties[i] <- propagated_uncertainty
    }
  }
}
mass_data[,1]=mass_data[,3]
mass_data$X.1[is.na(mass_data$X.1)] <- mass_data$propagated_uncertainties[is.na(mass_data$X.1)]
mass_data=mass_data[,1:2]
# Initialize a new dataframe to store the result
num_rows_md <- nrow(mass_data)
new_mass_data <- data.frame(matrix(NA, nrow = num_rows_md * 2, ncol = ncol(mass_data)))
colnames(new_mass_data) <- colnames(mass_data)

# Fill the new dataframe with data and empty rows
for (i in 1:num_rows_md) {
  new_mass_data[(i * 2 - 1), ] <- mass_data[i, ]  # Fill the data row
  # The even rows (i * 2) will remain NA by default
}

# Remove the last row if it is NA
if (all(is.na(new_mass_data[nrow(new_mass_data), ]))) {
  new_mass_data <- new_mass_data[-nrow(new_mass_data), ]
}

# Set the row names to run from 1 to the length of the data
rownames(new_mass_data) <- 1:nrow(new_mass_data)

# Assign the new dataframe back to mass_data
mass_data <- new_mass_data
# Insert an NA row at the beginning if the first element is not NA
if (!is.na(mass_data[1, 1])) {
  empty_row <- data.frame(X = NA, X.1 = NA)  # Create a new data frame row
  mass_data <- rbind(empty_row, mass_data)  # Bind the new row to the top
}
# Reset row names to ensure they are sequential
rownames(mass_data) <- 1:nrow(mass_data)
mass_data=as.matrix(mass_data)

mass_table[1:(nrow(mass_table)-1),4]=mass_data[,1]
mass_table[1:(nrow(mass_table)-1),5]=mass_data[,2]

# aerial mass calculations
mass_table[,6]=mass_table[,3]*mass_table[,4]
mass_table[,7]=mass_table[,5]*interval_thickness
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
    mass_table[i,10]=(mass_table[i-1,8]+mass_table[i+1,8])/2
  }
}
# mass depth section uncertainty 
for (i in 1:nrow(mass_table)-1){
  if (i %% 2 == 0) {
    mass_table[i,11]= sqrt(sum(mass_table[i+1,9]^2, mass_table[i-1,9]^2))/2
  }
}
# create the concentrations table and begin preliminary population
concentrations_table=array(data=NA, dim=c(((core_length*2)/interval_thickness), 14), dimnames=list(NULL, c("Mid Depth zi (cm)","Total Pb-210 dmp/g","Total Pb-210 Uncertainty u(Pb-210) dmp/g","Supported Pb-210 dmp/g","Supported Pb-210 Uncertainty u(sPb-210) dmp/g","Excess Pb-210 Ci dmp/g","Excess Pb-210 Uncertainty u(Ci) dmp/g","Excess Pb-210 C(i) dmp/g","Excess Pb-210 Uncertainty u(C(i)) dmp/g","Inventory delta Ai dpm/cm^2","Inventory Uncertainty u(delta Ai) dpm/cm^2","Mass Depth mi (g/cm^2)","ln(Ci)","u(ln(Ci))")))
concentrations_table[,1]=mass_table[1:nrow(concentrations_table),2]
concentrations_table[,12]=mass_table[1:nrow(concentrations_table),10]
# load in Pb-210 data, assign Alpha or Gamma methods, and fill activity data
alpha_or_gamma=readline(prompt = "Do you wish to enter Alpha data (True or False) ")
if (alpha_or_gamma==TRUE){
  print("select Alpha data")
  new_column_names_pb210_data <- c("X", "X.1", "X.2", "X.3")
  # Read the Excel file
  file_path2=file.choose()
  pb210_data <- read_excel(file_path2, col_names = TRUE, col_types = "numeric")
  if (ncol(pb210_data)==2){
    new_column_names_pb210_data=new_column_names_pb210_data[1:2]
  }
  colnames(pb210_data) <- new_column_names_pb210_data
  pb210_data=as.data.frame(pb210_data)
  # Bn=A*(C/A)^(n/N) exponential interpolation of activities
  for (i in 1:nrow(pb210_data)) {
    # Check if the current uncertainty is missing
    if (is.na(pb210_data[i, 2])) {
      # Get the first activity from the previous row
      first_activity_tot=pb210_data[(i-1),1]
      first_uncer_tot=pb210_data[(i-1),2]
      if (ncol(pb210_data)==4){
        first_activity_sup=pb210_data[(i-1),3]
        first_uncer_sup=pb210_data[(i-1),4]
      }
      # Initialize a variable to track the end of the missing data
      j = i
      # Loop to fill in missing uncertainties
      length_missing_data=1
      while (j <= nrow(pb210_data) && is.na(pb210_data[j, 2])) {
        j = j + 1
        length_missing_data=length_missing_data+1
      }
      length_missing_data=length_missing_data-1
      segments=length_missing_data+1
      filler_vec=vector("numeric", length_missing_data)
      filler_vec_uncer=vector("numeric", length_missing_data)
      if (ncol(pb210_data)==4){
        filler_vec2=vector("numeric", length_missing_data)
        filler_vec_uncer2=vector("numeric", length_missing_data)
      }
      # Now j is the first non-missing uncertainty after the missing values
      if (j <= nrow(pb210_data)) {
        second_activity_tot = pb210_data[j, 1]
        second_uncer_tot = pb210_data[j, 2]
        if (ncol(pb210_data)==4){
          second_activity_sup=pb210_data[j, 3]
          second_uncer_sup=pb210_data[j, 4]
        }
        # Check if both activities are numeric
        if (is.numeric(first_activity_tot) && is.numeric(second_activity_tot)) {
          # Calculate the filled activity
          for (k in 1:length_missing_data){
            filler_vec[k]=first_activity_tot*(second_activity_tot/first_activity_tot)^(k/segments)
            filler_vec_uncer[k]=filler_vec[k]*sqrt(((((segments-k)*(second_activity_tot/first_activity_tot)^(k/segments))/segments)*first_uncer_tot)^2+(((k*(second_activity_tot/first_activity_tot)^((k/segments)-1))/segments)*second_uncer_tot)^2)
          }
          # Fill the missing uncertainties
          pb210_data[i:(j - 1), 1] = filler_vec
          pb210_data[i:(j - 1), 2] = filler_vec_uncer
        }
      }
      # Move the index to j to continue checking for further missing values
      i = j - 1
    }
  }
  # Initialize a new dataframe to store the result
  num_rows_zi <- nrow(pb210_data)
  new_pb210_data <- data.frame(matrix(NA, nrow = num_rows_zi * 2, ncol = ncol(pb210_data)))
  colnames(new_pb210_data) <- colnames(pb210_data)
  
  # Fill the new dataframe with data and empty rows
  for (i in 1:num_rows_zi) {
    new_pb210_data[(i * 2 - 1), ] <- pb210_data[i, ]  # Fill the data row
    # The even rows (i * 2) will remain NA by default
  }
  
  # Remove the last row if it is NA
  if (all(is.na(new_pb210_data[nrow(new_pb210_data), ]))) {
    new_pb210_data <- new_pb210_data[-nrow(new_pb210_data), ]
  }
  
  # Set the row names to run from 1 to the length of the data
  rownames(new_pb210_data) <- 1:nrow(new_pb210_data)
  
  # Assign the new dataframe back to pb210_data
  pb210_data <- new_pb210_data
  # Insert an NA row at the beginning if the first element is not NA
  if (!is.na(pb210_data[1, 1])) {
    if (ncol(pb210_data)==2){
      empty_row <- data.frame(X = NA, X.1 = NA)  # Create a new data frame row
      pb210_data <- rbind(empty_row, pb210_data)  # Bind the new row to the top 
    }else{
      empty_row <- data.frame(X = NA, X.1 = NA, X.2=NA, X.3=NA)  # Create a new data frame row
      pb210_data <- rbind(empty_row, pb210_data)  # Bind the new row to the top
    }
    
  }
  # Reset row names to ensure they are sequential
  rownames(pb210_data) <- 1:nrow(pb210_data)
  pb210_data=as.matrix(pb210_data)
  
  concentrations_table[1:(nrow(concentrations_table)),2]=pb210_data[,1]
  concentrations_table[1:(nrow(concentrations_table)),3]=pb210_data[,2]
  if(ncol(pb210_data)==4){
    concentrations_table[1:(nrow(concentrations_table)),4]=pb210_data[,3]
    concentrations_table[1:(nrow(concentrations_table)),5]=pb210_data[,4]
  }
  
  accept_background_determination=FALSE
  
  #start of background auto detection
  pb210_data_auto_bkrng=na.omit(pb210_data)
  pb210_data_auto_bkrng=as.data.frame(pb210_data_auto_bkrng)
  pb210_data_auto_bkrng$X.2=log(pb210_data_auto_bkrng[,1])
  # Perform change point analysis
  cpt_result <- cpt.mean(pb210_data_auto_bkrng$X.2, method = "PELT", penalty = "SIC", minseglen = 3)
  # Get the change points
  change_points_bkgrnd <- cpts(cpt_result)
  if (length(change_points_bkgrnd) > 0) {
    suggested_background_index <- rownames(pb210_data_auto_bkrng)[change_points_bkgrnd[length(change_points_bkgrnd)]]
    cat("Suggested background activity begins at index:", suggested_background_index, "\n")
    print("Pb-210 Background Activity (dmp/g)")
    print(mean(pb210_data_auto_bkrng[change_points_bkgrnd[length(change_points_bkgrnd)]:length(pb210_data_auto_bkrng[,1]), 1]))
    print("Pb-210 Background Activity Uncertainty (dmp/g)")
    print(sd(pb210_data_auto_bkrng[change_points_bkgrnd[length(change_points_bkgrnd)]:length(pb210_data_auto_bkrng[,1]), 2]))
  } else {
    cat("No change points detected.\n")
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
    background_table[,1]=na.omit(concentrations_table[,2])
    background_table[,2]=na.omit(concentrations_table[,1])
    background_table[,3]=na.omit(concentrations_table[,3])
    background_table=as.data.frame(background_table)
    background_table_plot=ggplot(background_table, aes(x = x, y = y)) +
      geom_errorbar(aes(xmin = x - x_uncertainty, xmax = x + x_uncertainty), width = 0.2) +
      geom_point() +  # Add points
      labs(title = "Total Pb-210 vs Mid Depth", x = "Total Pb-210 (dmp/g)", y ="Mid Depth (cm)") +  # Add labels
      scale_x_continuous(breaks = seq(min(background_table$x), max(background_table$x), by = breaks_by_value(background_table, "x")), labels=scales::number_format(accuracy=accuracy_finder(background_table[,3]))) +  # Set x-axis dynamically
      scale_y_reverse(breaks = seq(min(background_table$y), max(background_table$y), by = 4)) +  # Set y-axis dynamically
      theme_minimal()  # Use a minimal theme
    print(background_table_plot)
    first_background=as.numeric(readline(prompt = "At what index does background begin? "))
    background_value=mean(na.omit(concentrations_table[first_background:nrow(concentrations_table),2]))
    check_for_negs=which((concentrations_table[1:(first_background-1),2]-background_value)<0)[1]
    check_background_length=length(concentrations_table[first_background:nrow(concentrations_table)])
    if (is.na(check_for_negs)==FALSE){
      accept_new_background=FALSE
      while(accept_new_background==FALSE){
        cat(sprintf("poor background selection detected! negative excess Pb-210 activity at index: %f\n",check_for_negs))
        accept_new_background=readline(prompt = "Accept new background determination from this index? (True) ")
        background_value=mean(na.omit(concentrations_table[check_for_negs:nrow(concentrations_table),2]))
        background_value_uncer=sd(na.omit(concentrations_table[check_for_negs:nrow(concentrations_table),2]))
        first_background=check_for_negs
      }
    }else{
      background_value=mean(na.omit(concentrations_table[first_background:nrow(concentrations_table),2]))
      background_value_uncer=sd(na.omit(concentrations_table[first_background:nrow(concentrations_table),2]))
    }
    if (check_background_length<5){
      print("Background selection does not contain enough points! Make a new selection")
      accept_background_determination=FALSE
    }else{
    print("Background activity (dmp/g)")
    print(background_value)
    print("Background activity uncertainty (dmp/g)")
    print(background_value_uncer)
    accept_background_determination=readline(prompt = "Accept background determination? (True or False) ")
    }
  }
  for (i in 1:nrow(concentrations_table)){
    if (i %% 2 == 0) {
      concentrations_table[i,4]=background_value
      concentrations_table[i,5]=background_value_uncer
    }
  }
  concentrations_table=concentrations_table[1:first_background-1,]
  mass_table=mass_table[1:first_background-1,]
} else {
  print("select Gamma data")
  # Define the path to your Excel file
  new_column_names_pb210_data <- c("X", "X.1", "X.2", "X.3")
  # Read the Excel file
  file_path2=file.choose()
  pb210_data <- read_excel(file_path2, col_names = TRUE, col_types = "numeric")
  colnames(pb210_data) <- new_column_names_pb210_data
  pb210_data=as.data.frame(pb210_data)
  
  # Perform linear interpolation on the supported values
  pb210_data$interpolated_X.2 <- na.approx(pb210_data$X.2, na.rm = FALSE)
  
  # Initialize a vector for propagated uncertainties
  pb210_data$propagated_uncertainties_X.2 <- NA
  
  # Loop through the interpolated supported values to calculate uncertainties
  for (i in 1:nrow(pb210_data)) {
    # Only calculate propagated uncertainty if the current uncertainty is NA
    if (is.na(pb210_data$X.3[i])) {
      # Find the nearest known X on both sides
      left_index <- i - 1
      while (left_index > 0 && (is.na(pb210_data$X.2[left_index]) || is.na(pb210_data$X.3[left_index]))) {
        left_index <- left_index - 1
      }
      
      right_index <- i + 1
      while (right_index < nrow(pb210_data) && (is.na(pb210_data$X.2[right_index]) || is.na(pb210_data$X.3[right_index]))) {
        right_index <- right_index + 1
      }
      
      # Check if we found valid indices
      if (left_index > 0 && right_index < nrow(pb210_data)) {
        delta_y1 <- pb210_data$X.3[left_index]  # Use the correct column for uncertainties
        x1 <- left_index
        delta_y2 <- pb210_data$X.3[right_index]  # Use the correct column for uncertainties
        x2 <- right_index
        partial_der=(i-x1)/(x2-x1)
        
        # Calculate propagated uncertainty
        propagated_uncertainty <- sqrt(((1-partial_der)*delta_y1)^2 + (partial_der * delta_y2)^2)
        pb210_data$propagated_uncertainties_X.2[i] <- propagated_uncertainty
      }
    }
  }
  pb210_data[,3]=pb210_data[,5]
  pb210_data$X.3[is.na(pb210_data$X.3)] <- pb210_data$propagated_uncertainties[is.na(pb210_data$X.3)]
  pb210_data=pb210_data[,1:4]
  # Bn=A*(C/A)^(n/N) exponential interpolation of activities
  for (i in 1:nrow(pb210_data)) {
    # Check if the current uncertainty is missing
    if (is.na(pb210_data[i, 2])) {
      # Get the first activity from the previous row
      first_activity_tot=pb210_data[(i-1),1]
      first_uncer_tot=pb210_data[(i-1),2]
      # Initialize a variable to track the end of the missing data
      j = i
      # Loop to fill in missing uncertainties
      length_missing_data=1
      while (j <= nrow(pb210_data) && is.na(pb210_data[j, 2])) {
        j = j + 1
        length_missing_data=length_missing_data+1
      }
      length_missing_data=length_missing_data-1
      segments=length_missing_data+1
      filler_vec=vector("numeric", length_missing_data)
      filler_vec_uncer=vector("numeric", length_missing_data)
      if (ncol(pb210_data)==4){
        filler_vec2=vector("numeric", length_missing_data)
        filler_vec_uncer2=vector("numeric", length_missing_data)
      }
      # Now j is the first non-missing uncertainty after the missing values
      if (j <= nrow(pb210_data)) {
        second_activity_tot = pb210_data[j, 1]
        second_uncer_tot = pb210_data[j, 2]
        # Check if both activities are numeric
        if (is.numeric(first_activity_tot) && is.numeric(second_activity_tot)) {
          # Calculate the filled activity
          for (k in 1:length_missing_data){
            filler_vec[k]=first_activity_tot*(second_activity_tot/first_activity_tot)^(k/segments)
            filler_vec_uncer[k]=filler_vec[k]*sqrt(((((segments-k)*(second_activity_tot/first_activity_tot)^(k/segments))/segments)*first_uncer_tot)^2+(((k*(second_activity_tot/first_activity_tot)^((k/segments)-1))/segments)*second_uncer_tot)^2)
          }
          # Fill the missing uncertainties
          pb210_data[i:(j - 1), 1] = filler_vec
          pb210_data[i:(j - 1), 2] = filler_vec_uncer
        }
      }
      # Move the index to j to continue checking for further missing values
      i = j - 1
    }
  }
  # Initialize a new dataframe to store the result
  num_rows_zi <- nrow(pb210_data)
  new_pb210_data <- data.frame(matrix(NA, nrow = num_rows_zi * 2, ncol = ncol(pb210_data)))
  colnames(new_pb210_data) <- colnames(pb210_data)
  
  # Fill the new dataframe with data and empty rows
  for (i in 1:num_rows_zi) {
    new_pb210_data[(i * 2 - 1), ] <- pb210_data[i, ]  # Fill the data row
    # The even rows (i * 2) will remain NA by default
  }
  
  # Remove the last row if it is NA
  if (all(is.na(new_pb210_data[nrow(new_pb210_data), ]))) {
    new_pb210_data <- new_pb210_data[-nrow(new_pb210_data), ]
  }
  
  # Set the row names to run from 1 to the length of the data
  rownames(new_pb210_data) <- 1:nrow(new_pb210_data)
  
  # Assign the new dataframe back to pb210_data
  pb210_data <- new_pb210_data
  # Insert an NA row at the beginning if the first element is not NA
  if (!is.na(pb210_data[1, 1])) {
    if (ncol(pb210_data)==2){
      empty_row <- data.frame(X = NA, X.1 = NA)  # Create a new data frame row
      pb210_data <- rbind(empty_row, pb210_data)  # Bind the new row to the top 
    }else{
      empty_row <- data.frame(X = NA, X.1 = NA, X.2=NA, X.3=NA)  # Create a new data frame row
      pb210_data <- rbind(empty_row, pb210_data)  # Bind the new row to the top
    }
    
  }
  # Reset row names to ensure they are sequential
  rownames(pb210_data) <- 1:nrow(pb210_data)
  pb210_data=as.matrix(pb210_data)
  
  concentrations_table[1:(nrow(concentrations_table)),2]=pb210_data[,1]
  concentrations_table[1:(nrow(concentrations_table)),3]=pb210_data[,2]
  if(ncol(pb210_data)==4){
    concentrations_table[1:(nrow(concentrations_table)),4]=pb210_data[,3]
    concentrations_table[1:(nrow(concentrations_table)),5]=pb210_data[,4]
  }
}
concentrations_table[,6]=concentrations_table[,2]-concentrations_table[,4]
negative_index=which(concentrations_table[,6]<0)
if (length(negative_index) > 0) {
  concentrations_table=concentrations_table[1:(negative_index[1] - 1),]
  mass_table=mass_table[1:(negative_index[1] - 1),]
}
if (is.na(concentrations_table[nrow(concentrations_table),1])==TRUE){
  concentrations_table=concentrations_table[1:(nrow(concentrations_table)-1),1:ncol(concentrations_table)]
}
# excess pb-210 uncertainty
for (i in 1:nrow(concentrations_table)){
  if (i %% 2 == 0) {
    concentrations_table[i,7]=sqrt(sum(concentrations_table[i,3]^2, concentrations_table[i,5]^2))
  }
}
# inventory, inventory uncertainty, mass depth, ln(Ci), u(ln(Ci)) calculations

for (i in 1:nrow(concentrations_table)){
  if (i %% 2 == 0) {
    concentrations_table[i,10]=concentrations_table[i,6]*mass_table[i,6]
    concentrations_table[i,11]=concentrations_table[i,10]*sqrt(sum((concentrations_table[i,7]/concentrations_table[i,6])^2,(mass_table[i,7]/mass_table[i,6])^2))
    concentrations_table[i,13]=log(concentrations_table[i,6])
    concentrations_table[i,14]=concentrations_table[i,7]/concentrations_table[i,6]
  }
}

#COCA calculations
C0CA_table=array(data=NA, dim=c(10,3), dimnames=list(NULL, c("zi","mi","ln(Ci)")))
C0CA_table[,1]=na.omit(mass_table[1:20,2])
C0CA_table[,2]=na.omit(concentrations_table[1:20,12])
C0CA_table[,3]=na.omit(concentrations_table[1:20,13])

C0CA_model_list=list()
for (i in 3:10) {
  # Define the extent of the data to use for the model
  extent <- i  # This will give 3, 4, ..., 10
  
  # Extract the x and y values from the array
  x_values <- C0CA_table[1:extent, 2]  # Column 2 for x values
  y_values <- C0CA_table[1:extent, 3]  # Column 3 for y values
  
  # Fit the linear model and store it in the list
  C0CA_model_list[[i-2]] <- lm(y_values ~ x_values)
}
# Step 1: Extract R^2 values
r_squared_values <- sapply(C0CA_model_list, function(model) summary(model)$adj.r.squared)

# Step 2: Identify the index of the model with the highest R^2
best_model_index <- which.max(r_squared_values)

# Step 3: Retrieve the corresponding model
best_model <- C0CA_model_list[[best_model_index]]

# C(i)
for (i in 1:nrow(concentrations_table)-1){
  if (i==1){
    concentrations_table[i,8]=exp(best_model$coefficients[1])
  }else{
    if (i %% 2 == 1) {
      concentrations_table[i,8]=(concentrations_table[i+1,6]+concentrations_table[i-1,6])/2
    }
  }
}
# C(i) uncertainty (hardcoded first value may require adjustment)
for (i in 1:nrow(concentrations_table)-1){
  if (i==1){
    concentrations_table[i,9]=concentrations_table[i,8]*(summary(best_model)$coefficients[1, 2])
  }else{
    if (i %% 2 == 1) {
      concentrations_table[i,9]= sqrt(sum(concentrations_table[i+1,7]^2, concentrations_table[i-1,7]^2))/2
    }
  }
}
#save off copies of background-corrected mass and concentration data, accompanying plots
# find your destination folder
print("select destination folder")
destination_folder <- choose.dir(default = "", caption = "Select Destination Folder")

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
folder_name = "Background_corrected_mass_concentration_data"

# Create the folder in the working directory
dir.create(folder_name)

# Define the file path where you want to save the file
mass_file_path <- file.path(folder_name, "mass_table.csv")

# Save the data frame as a CSV file in the new folder
write.csv(mass_table, mass_file_path, row.names = FALSE, na = "")

# Define the file path where you want to save the file
concentrations_file_path <- file.path(folder_name, "concentrations_data.csv")

# Save the data frame as a CSV file in the new folder
write.csv(concentrations_table, concentrations_file_path, row.names = FALSE, na = "")

#C0CA save
model_summaries <- data.frame(Model = character(),
                              Coefficient = character(),
                              Estimate = numeric(),
                              Std.Error = numeric(),
                              t.value = numeric(),
                              Pr.t = numeric(),
                              Adj.R.squared = numeric(),
                              stringsAsFactors = FALSE)

# Extract information from each model and store it in the data frame
for (model_index in seq_along(C0CA_model_list)) {
  model <- C0CA_model_list[[model_index]]
  
  # Check if the model is valid
  if (!is.null(model)) {
    summary_model <- summary(model)
    
    # Extract coefficients
    coefficients <- summary_model$coefficients
    for (i in 1:nrow(coefficients)) {
      model_summaries <- rbind(model_summaries, 
                               data.frame(Model = paste("Model", model_index, "Indices 1-",model_index+2),  # Adjusting index for naming
                                          Coefficient = rownames(coefficients)[i],
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
write_xlsx(model_summaries, path = "Background_corrected_mass_concentration_data/C0CA_model_summaries.xlsx")

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
y_breaks <- seq(min(dbd_plotting_table$y), max(dbd_plotting_table$y), by = breaks_by_value(dbd_plotting_table, "y"))
# Create the plot
dbd_plot = ggplot(dbd_plotting_table, aes(x = x, y = y)) +
  geom_point() +  # Add points
  geom_errorbar(aes(ymin = y - y_uncertainty, ymax = y + y_uncertainty), width = 0.2) +  # Add error bars
  labs(title = "Dry Bulk Density vs Depth", x = "Mid Depth (cm)", y = "Dry Bulk Density (g/cm^3)") +  # Add labels
  scale_x_continuous(breaks = seq(min(dbd_plotting_table$x), max(dbd_plotting_table$x), by = interval_thickness*2), labels=scales::number_format(accuracy=0.5)) +  # Set x-axis dynamically
  scale_y_continuous(breaks = y_breaks,  # Set y-axis breaks
                     labels = scales::number_format(accuracy = accuracy_finder(dbd_plotting_table[, 3]))) +  # Set y-axis labels
  theme_minimal()  # Use a minimal theme
#by = breaks_by_value(tpb210_plotting_table, "x")), labels=scales::number_format(accuracy=accuracy_finder(tpb210_plotting_table[,3]))) +  # Set x-axis dynamically
print(dbd_plot)
dev.off()
# excess pb-210
pdf(file.path(folder_name, "excess_pb210.pdf"), width = 8, height = 6)
epb210_plotting_table=array(data=NA, dim=c(nrow(mass_table)/2,3), dimnames=list(NULL, c("x","y","x_uncertainty")))
epb210_plotting_table[,1]=na.omit(concentrations_table[,6])
epb210_plotting_table[,2]=na.omit(concentrations_table[,1])
epb210_plotting_table[,3]=na.omit(concentrations_table[,7])
epb210_plotting_table=as.data.frame(epb210_plotting_table)
epb210_plot=ggplot(epb210_plotting_table, aes(x = x, y = y)) +
  geom_point() +  # Add points
  geom_errorbar(aes(xmin = x - x_uncertainty, xmax = x + x_uncertainty), width = 0.2) +  # Add error bars
  labs(title = "Excess Pb-210 vs Depth", x = "Excess Pb-210 (dmp/g)", y ="Mid Depth (cm)") +  # Add labels
  scale_x_continuous(breaks = seq(min(epb210_plotting_table$x), max(epb210_plotting_table$x), by = breaks_by_value(epb210_plotting_table, "x")), labels=scales::number_format(accuracy=accuracy_finder(epb210_plotting_table[,3]))) +  # Set x-axis dynamically
  scale_y_reverse(breaks = seq(min(epb210_plotting_table$y), max(epb210_plotting_table$y), by = interval_thickness*2), labels=scales::number_format(accuracy=0.5)) +  # Set y-axis dynamically
  theme_minimal()  # Use a minimal theme
print(epb210_plot)
dev.off()
# supported pb-210
pdf(file.path(folder_name, "supported_pb210.pdf"), width = 8, height = 6)
spb210_plotting_table=array(data=NA, dim=c(nrow(mass_table)/2,3), dimnames=list(NULL, c("x","y","x_uncertainty")))
spb210_plotting_table[,1]=na.omit(concentrations_table[,4])
spb210_plotting_table[,2]=na.omit(concentrations_table[,1])
spb210_plotting_table[,3]=na.omit(concentrations_table[,5])
spb210_plotting_table=as.data.frame(spb210_plotting_table)
spb210_plot = ggplot(spb210_plotting_table, aes(x = x, y = y)) +
  geom_point() +  # Add points
  geom_errorbar(aes(xmin = x - x_uncertainty, xmax = x + x_uncertainty), width = 0.2) +  # Add error bars
  labs(title = "Supported Pb-210 vs Depth", x = "Supported Pb-210 (dmp/g)", y = "Mid Depth (cm)") +  # Add labels
  scale_x_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 0.5),labels=scales::number_format(accuracy=accuracy_finder(epb210_plotting_table[,3]))) +  # Set x-axis limits and increments of 0.1
  scale_y_reverse(breaks = seq(min(spb210_plotting_table$y), max(spb210_plotting_table$y), by = interval_thickness*2), labels=scales::number_format(accuracy=0.5)) +  # Set x-axis dynamically
  theme_minimal()  # Use a minimal theme
print(spb210_plot)
dev.off()
# total pb-210
pdf(file.path(folder_name, "total_pb210.pdf"), width = 8, height = 6)
tpb210_plotting_table=array(data=NA, dim=c(nrow(mass_table)/2,3), dimnames=list(NULL, c("x","y","x_uncertainty")))
tpb210_plotting_table[,1]=na.omit(concentrations_table[,2])
tpb210_plotting_table[,2]=na.omit(concentrations_table[,1])
tpb210_plotting_table[,3]=na.omit(concentrations_table[,3])
tpb210_plotting_table=as.data.frame(tpb210_plotting_table)
tpb210_plot=ggplot(tpb210_plotting_table, aes(x = x, y = y)) +
  geom_point() +  # Add points
  geom_errorbar(aes(xmin = x - x_uncertainty, xmax = x + x_uncertainty), width = 0.2) +  # Add error bars
  labs(title = "Total Pb-210 vs Depth", x = "Total Pb-210 (dmp/g)", y ="Mid Depth (cm)") +  # Add labels
  scale_x_continuous(breaks = seq(min(tpb210_plotting_table$x), max(tpb210_plotting_table$x), by = breaks_by_value(tpb210_plotting_table, "x")), labels=scales::number_format(accuracy=accuracy_finder(tpb210_plotting_table[,3]))) +  # Set x-axis dynamically
  scale_y_reverse(breaks = seq(min(tpb210_plotting_table$y), max(tpb210_plotting_table$y), by = interval_thickness*2), labels=scales::number_format(accuracy=0.5)) +  # Set x-axis dynamically
  theme_minimal()  # Use a minimal theme
print(tpb210_plot)
dev.off()
# total log scale pb-210
pdf(file.path(folder_name, "excess_log_pb210.pdf"), width = 8, height = 6)
tlogpb210_plotting_table=array(data=NA, dim=c(nrow(mass_table)/2,3), dimnames=list(NULL, c("x","y","x_uncertainty")))
tlogpb210_plotting_table[,1]=na.omit(concentrations_table[,13])
tlogpb210_plotting_table[,2]=na.omit(concentrations_table[,1])
tlogpb210_plotting_table[,3]=na.omit(concentrations_table[,14])
tlogpb210_plotting_table=as.data.frame(tlogpb210_plotting_table)
tlogpb210_plot=ggplot(tlogpb210_plotting_table, aes(x = x, y = y)) +
  geom_point() +  # Add points
  geom_errorbar(aes(xmin = x - x_uncertainty, xmax = x + x_uncertainty), width = 0.2) +  # Add error bars
  labs(title = "Excess Log-scale Pb-210 vs Depth", x = "Excess Log-scale Pb-210 (dmp/g)", y ="Mid Depth (cm)") +  # Add labels
  scale_x_continuous(breaks = seq(min(tlogpb210_plotting_table$x), max(tlogpb210_plotting_table$x), by = breaks_by_value(tlogpb210_plotting_table, "x")), labels=scales::number_format(accuracy=accuracy_finder(tlogpb210_plotting_table[,3]))) +  # Set x-axis dynamically
  scale_y_reverse(breaks = seq(min(tlogpb210_plotting_table$y), max(tlogpb210_plotting_table$y), by = interval_thickness*2), labels=scales::number_format(accuracy=0.5)) +  # Set x-axis dynamically
  theme_minimal()  # Use a minimal theme
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
epb210_md_plot=ggplot(epb210_md_plotting_table, aes(x = x, y = y)) +
  geom_point() +  # Add points
  geom_errorbar(aes(xmin = x - x_uncertainty, xmax = x + x_uncertainty), width = 0.2) +  # Add error bars
  labs(title = "Excess Pb-210 vs Cumulative Mass Depth", x = "Excess Pb-210 (dmp/g)", y ="Cumulative Mass Depth (g/cm^2)") +  # Add labels
  scale_x_continuous(breaks = seq(min(epb210_md_plotting_table$x), max(epb210_md_plotting_table$x), by = breaks_by_value(epb210_md_plotting_table, "x")), labels=scales::number_format(accuracy=accuracy_finder(epb210_md_plotting_table[,3]))) +  # Set x-axis dynamically
  scale_y_reverse(breaks = seq(min(epb210_md_plotting_table$y), max(epb210_md_plotting_table$y), by = breaks_by_value(mass_depth_table, "Mass Depth mi (g/cm^2)")), labels=scales::number_format(accuracy=accuracy_finder(mass_depth_table[,2]))) +  # Set y-axis dynamically
  theme_minimal()  # Use a minimal theme
print(epb210_md_plot)
dev.off()
# supported pb-210 mass depth
pdf(file.path(folder_name, "supported_pb210_md.pdf"), width = 8, height = 6)
spb210_md_plotting_table=array(data=NA, dim=c(nrow(mass_table)/2,3), dimnames=list(NULL, c("x","y","x_uncertainty")))
spb210_md_plotting_table[,1]=na.omit(concentrations_table[,4])
spb210_md_plotting_table[,2]=na.omit(concentrations_table[,12])
spb210_md_plotting_table[,3]=na.omit(concentrations_table[,5])
spb210_md_plotting_table=as.data.frame(spb210_md_plotting_table)
spb210_md_plot = ggplot(spb210_md_plotting_table, aes(x = x, y = y)) +
  geom_point() +  # Add points
  geom_errorbar(aes(xmin = x - x_uncertainty, xmax = x + x_uncertainty), width = 0.2) +  # Add error bars
  labs(title = "Supported Pb-210 vs Cumulative Mass Depth", x = "Supported Pb-210 (dmp/g)", y = "Cumulative Mass Depth (g/cm^2)") +  # Add labels
  scale_x_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 0.5),labels=scales::number_format(accuracy=accuracy_finder(epb210_plotting_table[,3]))) +
  scale_y_reverse(breaks = seq(min(spb210_md_plotting_table$y), max(spb210_md_plotting_table$y), by = breaks_by_value(mass_depth_table, "Mass Depth mi (g/cm^2)")), labels=scales::number_format(accuracy=accuracy_finder(mass_depth_table[,2]))) +  # Set y-axis dynamically
  theme_minimal()  # Use a minimal theme
print(spb210_md_plot)
dev.off()
# total pb-210 mass depth 
pdf(file.path(folder_name, "total_pb210_md.pdf"), width = 8, height = 6)
tpb210_md_plotting_table=array(data=NA, dim=c(nrow(mass_table)/2,3), dimnames=list(NULL, c("x","y","x_uncertainty")))
tpb210_md_plotting_table[,1]=na.omit(concentrations_table[,2])
tpb210_md_plotting_table[,2]=na.omit(concentrations_table[,12])
tpb210_md_plotting_table[,3]=na.omit(concentrations_table[,3])
tpb210_md_plotting_table=as.data.frame(tpb210_md_plotting_table)
tpb210_md_plot=ggplot(tpb210_md_plotting_table, aes(x = x, y = y)) +
  geom_point() +  # Add points
  geom_errorbar(aes(xmin = x - x_uncertainty, xmax = x + x_uncertainty), width = 0.2) +  # Add error bars
  labs(title = "Total Pb-210 vs Cumulative Mass Depth", x = "Total Pb-210 (dmp/g)", y ="Cumulative Mass Depth (g/cm^2)") +  # Add labels
  scale_x_continuous(breaks = seq(min(tpb210_md_plotting_table$x), max(tpb210_md_plotting_table$x), by = breaks_by_value(tpb210_md_plotting_table, "x")), labels=scales::number_format(accuracy=accuracy_finder(tpb210_md_plotting_table[,3]))) +  # Set x-axis dynamically
  scale_y_reverse(breaks = seq(min(tpb210_md_plotting_table$y), max(tpb210_md_plotting_table$y), by = breaks_by_value(mass_depth_table, "Mass Depth mi (g/cm^2)")), labels=scales::number_format(accuracy=accuracy_finder(mass_depth_table[,2]))) +  # Set y-axis dynamically
  theme_minimal()  # Use a minimal theme
print(tpb210_md_plot)
dev.off()
# Excess log scale pb-210 mass depth
pdf(file.path(folder_name, "excess_log_pb210_md.pdf"), width = 8, height = 6)
tlogpb210_md_plotting_table=array(data=NA, dim=c(nrow(mass_table)/2,3), dimnames=list(NULL, c("x","y","x_uncertainty")))
tlogpb210_md_plotting_table[,1]=na.omit(concentrations_table[,13])
tlogpb210_md_plotting_table[,2]=na.omit(concentrations_table[,12])
tlogpb210_md_plotting_table[,3]=na.omit(concentrations_table[,14])
tlogpb210_md_plotting_table=as.data.frame(tlogpb210_md_plotting_table)
tlogpb210_md_plot=ggplot(tlogpb210_md_plotting_table, aes(x = x, y = y)) +
  geom_point() +  # Add points
  geom_errorbar(aes(xmin = x - x_uncertainty, xmax = x + x_uncertainty), width = 0.2) +  # Add error bars
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Excess Log-scale Pb-210 vs Cumulative Mass Depth", x = "Excess Log-scale Pb-210 (dmp/g)", y ="Cumulative Mass Depth (g/cm^2)") +  # Add labels
  scale_x_continuous(breaks = seq(min(tlogpb210_md_plotting_table$x), max(tlogpb210_md_plotting_table$x), by = breaks_by_value(tlogpb210_md_plotting_table, "x")), labels=scales::number_format(accuracy=accuracy_finder(tlogpb210_md_plotting_table[,3]))) +  # Set x-axis dynamically
  scale_y_reverse(breaks = seq(min(tlogpb210_md_plotting_table$y), max(tlogpb210_md_plotting_table$y), by = breaks_by_value(mass_depth_table, "Mass Depth mi (g/cm^2)")), labels=scales::number_format(accuracy=accuracy_finder(mass_depth_table[,2]))) +  # Set y-axis dynamically
  theme_minimal()  # Use a minimal theme
print(tlogpb210_md_plot)
dev.off()
#make constants table to save
constants_table=array(data =NA, dim=c(14, 2), dimnames=list(NULL,c("constant","value")))
constants_table[,1]=c("day","month","year","full numeric year","core diameter (cm)","core diameter uncertainty (cm)","surface area (cm^2)","surface area uncertainty (cm^2)","interval thickness (cm)","full core length (cm)","Pb-210 half-life (yr)","Pb-210 half-life uncertainty (yr)","Pb-210 disintegration constant (yr^-1)","Pb-210 disinetragtion constant uncertainty (yr^-1)")
constants_table[1,2]=day
constants_table[2,2]=month
constants_table[3,2]=year
constants_table[4,2]=numeric_date
constants_table[5,2]=internal_diameter
constants_table[6,2]=internal_diameter_uncer
constants_table[7,2]=surface_area
constants_table[8,2]=surface_area_uncer
constants_table[9,2]=interval_thickness
constants_table[10,2]=core_length
constants_table[11,2]=half_life
constants_table[12,2]=half_life_uncer
constants_table[13,2]=decay_const
constants_table[14,2]=decay_const_uncer
constants_table=as.data.frame(constants_table)
constants_table$value=as.numeric(constants_table$value)
constants_file_path <- file.path(folder_name, "constants_table.csv")
write.csv(constants_table, constants_file_path, row.names=FALSE)
#start of CFCS section
CFCS_table=array(data=NA, dim=c(nrow(concentrations_table)/2,6), dimnames=list(NULL,c("zi","mi","ln(Ci)","Time (yr)","Date (yr)", "Time Uncertainty (yr)")))
CFCS_table[,1]=na.omit(concentrations_table[,1])
CFCS_table[,2]=na.omit(concentrations_table[,12])
CFCS_table[,3]=na.omit(concentrations_table[,13])
CFCS_table=as.data.frame(CFCS_table)
CFCS_modeling_complete=FALSE
SAZ_depth=0
  surface_table=array(data=NA, dim=c(nrow(mass_table)/2,3), dimnames=list(NULL, c("x","y", "x_uncertainty")))
  surface_table[,1]=na.omit(concentrations_table[,6])
  surface_table[,2]=na.omit(concentrations_table[,1])
  surface_table[,3]=na.omit(concentrations_table[,7])
  surface_table=as.data.frame(surface_table)
  surface_table_plot=ggplot(surface_table, aes(x = x, y = y)) +
    geom_errorbar(aes(xmin = x - x_uncertainty, xmax = x + x_uncertainty), width = 0.2) +
    geom_point() +  # Add points
    labs(title = "Excess Pb-210 vs Mid Depth", x = "Excess Pb-210 (dmp/g)", y ="Mid Depth (cm)") +  # Add labels
    scale_x_continuous(breaks = seq(min(surface_table$x), max(surface_table$x), by = breaks_by_value(surface_table, "x")), labels=scales::number_format(accuracy=accuracy_finder(surface_table[,3]))) +  # Set x-axis dynamically
    scale_y_reverse(breaks = seq(min(surface_table$y), max(surface_table$y), by = 4)) +  # Set y-axis dynamically
    theme_minimal()  # Use a minimal theme
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
    colnames(temp_surface_table) = c("Excess Pb-210 (dmp/g)", "Mid Depth (cm)")
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
print("SAR modeling")
while (CFCS_modeling_complete==FALSE){
  SAR_model_plot=ggplot(CFCS_table, aes(x = CFCS_table[, 1], y = CFCS_table[, 3])) +
    geom_point(color = "blue", size = 2) +  # Customize point color and size
    labs(x = "Mid Depth (cm)", y = "log-scale Excess Pb-210 Activity (dpm/g)", title = "SAR Modeling") +
    scale_x_continuous(breaks = seq(0, max(CFCS_table[, 1]), by = interval_thickness*2), labels=scales::number_format(accuracy=0.5)) +  #fix this
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))  # Center the title
  print(SAR_model_plot)
  changepoint_detection_model = lm(`ln(Ci)` ~ zi, data = CFCS_table)
  number_of_models=as.numeric(readline(prompt = "How many models do you wish to try? "))
  number_of_changepoints=number_of_models-1
  changepoint_list=c()
  if (number_of_changepoints>0){
    for (i in 1:number_of_changepoints){
      changepoint_list[i]=as.numeric(readline(prompt = "Enter change-point "))
    }
    segmented_model <- segmented(changepoint_detection_model, seg.Z = ~ zi, psi = list(zi = changepoint_list))
    estimated_breakpoints <- summary(segmented_model)$psi[, "Est."]
    estimated_breakpoints=unlist(estimated_breakpoints, use.names = FALSE)
    SAR_indexes <- numeric(length(estimated_breakpoints))
    for (i in 1:length(estimated_breakpoints)) {
      
      SAR_indexes[i]<- which.min(abs(CFCS_table[,1]-estimated_breakpoints[i]))
    }
    print("estimated change point(s)")
    print(estimated_breakpoints)
    extents=c(1)
    extents <- c(extents, lapply(SAR_indexes, function(x) x))
    extents[length(extents)+1]=nrow(CFCS_table)
    print("model data values")
    model_data_vals=c()
    for (i in 1:length(extents)){
      model_data_vals[[i]]=CFCS_table[extents[[i]],1]
    }
    cat(unlist(model_data_vals), sep = " ", fill=TRUE)
    print("model data indexes")
    list_of_models=c()
    extents <- unlist(extents, use.names = FALSE)
    print(extents)
    for (i in 1:number_of_models){
      list_of_models[[i]]=lm(CFCS_table[extents[i]:extents[i+1],3]~CFCS_table[extents[i]:extents[i+1],1])
      print(paste("Model Adj. R square",i))
      print(summary(list_of_models[[i]])$adj.r.square)
      print(paste("Model SAR rate (cm/yr)",i))
      print(-decay_const/summary(list_of_models[[i]])$coefficients[2, 1])
      SAR_model_plot + geom_line(aes(y = fitted(list_of_models[[i]])), color = i, linewidth = 1)
    }
    lines_df <- data.frame()  # To hold all segments
    
    for (i in seq_along(list_of_models)) {
      model <- list_of_models[[i]]
      model_data <- model$model
      
      x_vals <- model_data[[2]]  # zi
      y_vals <- fitted(model)    # fitted ln(Ci)
      
      segment_df <- data.frame(x = x_vals, y = y_vals, group = as.factor(i))
      lines_df <- rbind(lines_df, segment_df)
    }
    print(ggplot(CFCS_table, aes(x = zi, y = `ln(Ci)`)) +
            geom_point(color = "blue") +
            labs(title = "SAR Modeling", x="Mid Depth (cm)", y="log-scale Excess Pb-210 Activity (dpm/g)") +
            geom_line(data = lines_df, aes(x = x, y = y, group = group, color = group), linewidth = 1.2, show.legend = TRUE) +
            labs(color = "Model"))
  }else{
    print("single model fit")
    print(summary(changepoint_detection_model)$adj.r.squared)
    print(SAR_model_plot + geom_line(aes(y = fitted(changepoint_detection_model)), color = "red", linewidth = 1))
    print("single model SAR rate (cm/yr)")
    print(-decay_const/summary(changepoint_detection_model)$coefficients[2, 1])
  }
  CFCS_manual_eval=readline(prompt = "Do you wish to manually set CFCS models? (True or False) ")
  if (CFCS_manual_eval==TRUE){
    SAR_model_plot=ggplot(CFCS_table, aes(x = CFCS_table[, 1], y = CFCS_table[, 3])) +
      geom_point(color = "blue", size = 2) +  # Customize point color and size
      labs(x = "Mid Depth (cm)", y = "log-scale Excess Pb-210 Activity (dpm/g)", title = "SAR Modeling") +
      scale_x_continuous(breaks = seq(0, max(CFCS_table[, 1]), by = interval_thickness*2), labels=scales::number_format(accuracy=0.5)) +  #fix this
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))  # Center the title
    print(SAR_model_plot)
    changepoint_detection_model = lm(`ln(Ci)` ~ zi, data = CFCS_table)
    number_of_models=as.numeric(readline(prompt = "How many models do you wish to try? "))
    number_of_changepoints=number_of_models-1
    changepoint_list=c()
    if (number_of_changepoints>0){
      for (i in 1:number_of_changepoints){
        changepoint_list[i]=as.numeric(readline(prompt = "Enter change-point "))
      }
      
    }
    extents=c(1, changepoint_list, length(CFCS_table[,1]))
    print("model data indexes")
    print(extents)
    print("model data values (cm)")
    model_data_vals=c()
    for (i in 1:length(extents)){
      model_data_vals[[i]]=CFCS_table[extents[[i]],1]
    }
    print(model_data_vals)
    list_of_models=c()
    for (i in 1:number_of_models){
      list_of_models[[i]]=lm(CFCS_table[extents[i]:extents[i+1],3]~CFCS_table[extents[i]:extents[i+1],1])
      print(paste("Model Adj. R square",i))
      print(summary(list_of_models[[i]])$adj.r.square)
      print(paste("Model SAR rate (cm/yr)",i))
      print(-decay_const/summary(list_of_models[[i]])$coefficients[2, 1])
      SAR_model_plot + geom_line(aes(y = fitted(list_of_models[[i]])), color = i, linewidth = 1)
    }
    lines_df <- data.frame()  # To hold all segments
    
    for (i in seq_along(list_of_models)) {
      model <- list_of_models[[i]]
      model_data <- model$model
      
      x_vals <- model_data[[2]]  # zi
      y_vals <- fitted(model)    # fitted ln(Ci)
      
      segment_df <- data.frame(x = x_vals, y = y_vals, group = as.factor(i))
      lines_df <- rbind(lines_df, segment_df)
    }
    print(ggplot(CFCS_table, aes(x = zi, y = `ln(Ci)`)) +
            geom_point(color = "blue") +
            labs(title = "SAR Modeling", x="Mid Depth (cm)", y="log-scale Excess Pb-210 Activity (dpm/g)") +
            geom_line(data = lines_df, aes(x = x, y = y, group = group, color = group), linewidth = 1.2, show.legend = TRUE) +
            labs(color = "Model"))
  }
  CFCS_modeling_complete=readline(prompt = "Do you wish to accept current CFCS model(s)? (True or False) ")
}

print("MAR modeling")
MAR_model_plot=ggplot(CFCS_table, aes(x = CFCS_table[, 2], y = CFCS_table[, 3])) +
  geom_point(color = "blue", size = 2) +  # Customize point color and size
  labs(x = "Cumulative Mass Depth (g/cm^2)", y = "log-scale Excess Pb-210 Activity (dpm/g)", title = "MAR Modeling") +
  scale_x_continuous(breaks = seq(0, max(CFCS_table[, 2]), by = 1)) +  # Set x-axis breaks
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title
print(MAR_model_plot)
changepoint_detection_model_MAR = lm(`ln(Ci)` ~ mi, data = CFCS_table)

if (number_of_changepoints>0){
  list_of_models_MAR=c()
  print("model data values")
  model_data_vals_MAR=c()
  for (i in 1:length(extents)){
    model_data_vals_MAR[[i]]=CFCS_table[extents[[i]],2]
  }
  cat(unlist(model_data_vals_MAR), sep = " ", fill=TRUE)
  print("model data indexes")
  print(extents)
  for (i in 1:number_of_models){
    list_of_models_MAR[[i]]=lm(CFCS_table[extents[i]:extents[i+1],3]~CFCS_table[extents[i]:extents[i+1],2])
    print(paste("Model Adj. R square",i))
    print(summary(list_of_models_MAR[[i]])$adj.r.square)
    print(paste("Model MAR rate (g/(cm^2 yr)",i))
    print(-decay_const/summary(list_of_models_MAR[[i]])$coefficients[2, 1])
    MAR_model_plot + geom_line(aes(y = fitted(list_of_models_MAR[[i]])), color = i, linewidth = 1)
  }
  lines_df_MAR <- data.frame()  # To hold all segments
  
  for (i in seq_along(list_of_models_MAR)) {
    model_MAR <- list_of_models_MAR[[i]]
    model_data_MAR <- model_MAR$model
    
    x_vals_MAR <- model_data_MAR[[2]]  # mi
    y_vals_MAR <- fitted(model_MAR)    # fitted ln(Ci)
    
    segment_df_MAR <- data.frame(x = x_vals_MAR, y = y_vals_MAR, group = as.factor(i))
    lines_df_MAR <- rbind(lines_df_MAR, segment_df_MAR)
  }
  print(ggplot(CFCS_table, aes(x = mi, y = `ln(Ci)`)) +
          geom_point(color = "blue") +
          labs(title = "MAR Modeling", x="Cumulative Mass Depth (g/cm^2)", y="log-scale Pb-210 Activity (dpm/g)") +
          geom_line(data = lines_df_MAR, aes(x = x, y = y, group = group, color = group), linewidth = 1.2, show.legend = TRUE) +
          labs(color = "Model"))
  
}else{
  print("single model fit")
  print(summary(changepoint_detection_model_MAR)$adj.r.squared)
  print("single model MAR rate (g/(cm^2 yr)")
  print(-decay_const/summary(changepoint_detection_model_MAR)$coefficients[2, 1])
  print(MAR_model_plot + geom_line(aes(y = fitted(changepoint_detection_model)), color = "red", linewidth = 1))
  
}
#populate CFCS table
if (number_of_models>1){
  time_correction=0
  for (i in 1:length(list_of_models)){
    for (j in extents[i]:(extents[i+1])){
      if (i==1){
        CFCS_table[j,4]=(CFCS_table[j,1]-CFCS_table[extents[i],1])/(-decay_const/(summary(list_of_models[[i]])$coefficients[2, "Estimate"]))+time_correction
        CFCS_table[j,5]=numeric_date-CFCS_table[j,4]
      }else{
        CFCS_table[j,4]=(CFCS_table[j,1]-CFCS_table[(extents[i]-1),1])/(-decay_const/(summary(list_of_models[[i]])$coefficients[2, "Estimate"]))+time_correction
        CFCS_table[j,5]=numeric_date-CFCS_table[j,4]
      }
      if (j>=extents[i+1]-1 && j+1!=extents[i+1]){
        time_correction=CFCS_table[j-1,4]
      }
    }
  }
}else{
  for (i in 1:nrow(CFCS_table)){
    CFCS_table[i,4]=(CFCS_table[i,1]-CFCS_table[1,1])/(-decay_const/(summary(changepoint_detection_model)$coefficients[2, "Estimate"]))
    CFCS_table[i,5]=numeric_date-CFCS_table[i,4]
  }
}
# calculate the age uncertainty
if (number_of_models>=2){
for (i in 1:number_of_models){
  for (j in extents[[i]]:extents[[i+1]]){
    CFCS_table[j,6]=CFCS_table[j,4]*((-summary(list_of_models[[i]])$coefficients[2, 2])/(summary(list_of_models[[i]])$coefficients[2,1]))
  }
}
}else{ 
  for (i in 1:nrow(CFCS_table)){
    CFCS_table[i,6]=CFCS_table[i,4]*((-summary(changepoint_detection_model)$coefficients[2, 2])/(summary(changepoint_detection_model)$coefficients[2,1]))
  }
}

folder_name2 = "CFCS modeling data"
dir.create(folder_name2)
CFCS_table_path <- file.path(folder_name2, "CFCS_table.csv")
write.csv(CFCS_table, CFCS_table_path, row.names = FALSE, na = "")
# save the MAR/SAR rates and uncertainties to files
#MAR version
rownames_MAR_rates=c()
if (number_of_models==1){
MAR_rates = array(data=NA, dim=c(1, 2),dimnames=list(NULL, c("MAR (g/(cm^2 yr)", "uncertainty")))
MAR_rates[1,1]=-decay_const/summary(changepoint_detection_model_MAR)$coefficients[2, 1]
MAR_rates[1,2]=MAR_rates[1,1]*sqrt((-summary(changepoint_detection_model_MAR)$coefficients[2, 2]/summary(changepoint_detection_model_MAR)$coefficients[2, 1])^2+(decay_const_uncer/decay_const)^2)
rownames(MAR_rates)=c("single model")
}else{
MAR_rates = array(data=NA, dim=c(length(list_of_models), 2),dimnames=list(NULL, c("MAR (g/(cm^2 yr)", "uncertainty")))
for ( i in 1:length(list_of_models)){
  MAR_rates[i,1]=-decay_const/summary(list_of_models_MAR[[i]])$coefficients[2, 1]
  MAR_rates[i,2]=MAR_rates[i,1]*sqrt((-summary(list_of_models_MAR[[i]])$coefficients[2, 2]/summary(list_of_models_MAR[[i]])$coefficients[2, 1])^2+(decay_const_uncer/decay_const)^2)
  rownames_MAR_rates[i]=paste("Model",i)
}
rownames(MAR_rates)=rownames_MAR_rates
}
MAR_rates_path <- file.path(folder_name2, "MAR_rates.csv")
if (number_of_models>1){
  MAR_final_slope=summary(list_of_models_MAR[[length(list_of_models_MAR)]])$coefficients[2, "Estimate"]
  MAR_final_slope_err=summary(list_of_models_MAR[[length(list_of_models_MAR)]])$coefficients[2, 2]
}else{
  MAR_final_slope=summary(changepoint_detection_model_MAR)$coefficients[2, "Estimate"]
  MAR_final_slope_err=summary(changepoint_detection_model_MAR)$coefficients[2, 2]
}

write.csv(MAR_rates, MAR_rates_path, row.names = TRUE)
missing_inventory=(concentrations_table[(nrow(concentrations_table)),6]*(-decay_const/MAR_final_slope))/decay_const
MAR_final=(-decay_const/MAR_final_slope)
MAR_final_err=(-decay_const/MAR_final_slope)*sqrt((MAR_final_slope_err/MAR_final_slope)^2+(decay_const_uncer/decay_const)^2)
#SAR version
rownames_SAR_rates=c()
if (number_of_models==1){
  SAR_rates = array(data=NA, dim=c(1, 2),dimnames=list(NULL, c("SAR (cm/yr)", "uncertainty")))
  SAR_rates[1,1]=-decay_const/summary(changepoint_detection_model)$coefficients[2, 1]
  SAR_rates[1,2]=SAR_rates[1,1]*sqrt((-summary(changepoint_detection_model)$coefficients[2, 2]/summary(changepoint_detection_model)$coefficients[2, 1])^2+(decay_const_uncer/decay_const)^2)
  rownames(SAR_rates)=c("single model")
}else{
  SAR_rates = array(data=NA, dim=c(length(list_of_models), 2),dimnames=list(NULL, c("SAR (cm/yr)", "uncertainty")))
  for ( i in 1:length(list_of_models)){
    SAR_rates[i,1]=-decay_const/summary(list_of_models[[i]])$coefficients[2, 1]
    SAR_rates[i,2]=SAR_rates[i,1]*sqrt((-summary(list_of_models[[i]])$coefficients[2, 2]/summary(list_of_models[[i]])$coefficients[2, 1])^2+(decay_const_uncer/decay_const)^2)
    rownames_SAR_rates[i]=paste("Model",i)
  }
  rownames(SAR_rates)=rownames_SAR_rates
}
SAR_rates_path <- file.path(folder_name2, "SAR_rates.csv")
write.csv(SAR_rates, SAR_rates_path, row.names = TRUE)

##save changepoint models, if they are used
if (number_of_models==1){
  SAR_model_summary <- data.frame(Model = character(),
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
      SAR_model_summary <- rbind(SAR_model_summary, 
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
    if (nrow(SAR_model_summary) > 0) {
      SAR_model_summary$Model <- "Changepoint Detection Model"
    }
  } else {
    cat("The changepoint_detection_model is NULL or invalid.\n")
  }
  
  # Write the summary to an Excel file
  write_xlsx(SAR_model_summary, path = "CFCS modeling data/SAR_model_summary.xlsx")
  
  #now for the MAR version
  
  
  MAR_model_summary <- data.frame(Model = character(),
                                  Coefficient = character(),
                                  Estimate = numeric(),
                                  Std.Error = numeric(),
                                  t.value = numeric(),
                                  Pr.t = numeric(),
                                  Adj.R.squared = numeric(),
                                  stringsAsFactors = FALSE)
  
  # Check if the model is valid
  if (!is.null(changepoint_detection_model_MAR)) {
    summary_model <- summary(changepoint_detection_model_MAR)
    
    # Extract coefficients
    coefficients <- summary_model$coefficients
    for (i in 1:nrow(coefficients)) {
      MAR_model_summary <- rbind(MAR_model_summary, 
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
    if (nrow(MAR_model_summary) > 0) {
      MAR_model_summary$Model <- "Changepoint Detection Model"
    }
  } else {
    cat("The changepoint_detection_model is NULL or invalid.\n")
  }
  
  # Write the summary to an Excel file
  write_xlsx(MAR_model_summary, path = "CFCS modeling data/MAR_model_summary.xlsx")
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
  for (model_index in seq_along(list_of_models_MAR)) {
    model <- list_of_models_MAR[[model_index]]
    
    # Check if the model is valid
    if (!is.null(model)) {
      summary_model <- summary(model)
      
      # Extract coefficients
      coefficients <- summary_model$coefficients
      for (i in 1:nrow(coefficients)) {
        MAR_model_summaries <- rbind(MAR_model_summaries, 
                                     data.frame(Model = paste("Model", model_index, "Indicies",extents[model_index] , "-", extents[model_index+1]),  # Adjusting index for naming
                                                Coefficient = rownames(coefficients)[i],
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
## MAR lm summary save
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
  for (model_index in seq_along(list_of_models)) {
    model <- list_of_models[[model_index]]
    
    # Check if the model is valid
    if (!is.null(model)) {
      summary_model <- summary(model)
      
      # Extract coefficients
      coefficients <- summary_model$coefficients
      for (i in 1:nrow(coefficients)) {
        SAR_model_summaries <- rbind(SAR_model_summaries, 
                                     data.frame(Model = paste("Model", model_index, "Indicies",extents[model_index] , "-", extents[model_index+1]),  # Adjusting index for naming
                                                Coefficient = rownames(coefficients)[i],
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

# plot saves
#SAR
pdf(file.path(folder_name2, "SAR_models.pdf"), width = 8, height = 6)
if (number_of_models>1){
  print(ggplot(CFCS_table, aes(x = zi, y = `ln(Ci)`)) +
          geom_point(color = "blue") +
          labs(title = "SAR Modeling", x="Mid Depth (cm)", y="log-scale Pb-210 Activity (dpm/g)") +
          geom_line(data = lines_df, aes(x = x, y = y, group = group, color = group), linewidth = 1.2, show.legend = TRUE) +
          labs(color = "Model"))
}else{
  #SAR
  
  print(SAR_model_plot + geom_line(aes(y = fitted(changepoint_detection_model)), color = "red", linewidth = 1))
}
dev.off()

#MAR
pdf(file.path(folder_name2, "MAR_models.pdf"), width = 8, height = 6)
if (number_of_models>1){
  print(ggplot(CFCS_table, aes(x =mi, y = `ln(Ci)`)) +
          geom_point(color = "blue") +
          labs(title = "MAR Modeling", x="Cumulative Mass Depth (g/cm^2)", y="log-scale Pb-210 Activity (dpm/g)") +
          geom_line(data = lines_df_MAR, aes(x = x, y = y, group = group, color = group), linewidth = 1.2, show.legend = TRUE) +
          labs(color = "Model"))
}else{
  
  print(MAR_model_plot + geom_line(aes(y = fitted(changepoint_detection_model_MAR)), color = "red", linewidth = 1))
}
dev.off()

#age-depth
pdf(file.path(folder_name2, "age_depth.pdf"), width = 8, height = 6)
print(ggplot(CFCS_table, aes(x = `Date (yr)`, y = `zi`)) +
        geom_point(color = "blue") +
        geom_errorbar(aes(xmin = `Date (yr)` - `Time Uncertainty (yr)`, 
                           xmax = `Date (yr)` + `Time Uncertainty (yr)`, 
                           y = `zi`,), 
                       height = 0.2, color = "black") +  # Adjust height and color as needed
        labs(title = "Age vs Depth Model") +
        scale_y_reverse() +
        scale_x_reverse())
dev.off()

# age-mass depth
pdf(file.path(folder_name2, "age_mass_depth.pdf"), width = 8, height = 6)
print(ggplot(CFCS_table, aes(x = `Date (yr)`, y = `mi`)) +
        geom_point(color = "blue") +
        geom_errorbar(aes(xmin = `Date (yr)` - `Time Uncertainty (yr)`, 
                           xmax = `Date (yr)` + `Time Uncertainty (yr)`, 
                           y = `mi`), 
                       height = 0.2, color = "black") +  # Adjust height and color as needed
        labs(title = "Age vs Mass Depth Model") +
        scale_y_reverse() +
        scale_x_reverse())
dev.off()

if (surface_active_zone_check==TRUE){
  print("CFCS model assumptions violated! This core has a surface active zone, only the relative ages of the sediment layers below the surface active layer could be calculated, not absolute ages")
  CFCS_table[,5]=NA
  }

finish_CFCS=readline(prompt = "CFCS modeling complete. Continue to CF? (True) ")
#start of CF section
#create CF table and begin populating it
CF_table=array(data=NA, dim=c(nrow(mass_table),17), dimnames=list(NULL,c("Top of Interval z(i) (cm)", "Mid Depth zi (cm)","Inventory delta Ai dpm/cm^2","Inventory Uncertainty u(delta Ai) dpm/cm^2","Inventory Below A(i) dpm/cm^2","Inventory Below Uncertainty u(A(i)) dpm/cm^2","Age t(i) yr","Age Uncertainty u(t(i)) yr","Year CE","Mass Accumulation Rate r(i) g/(cm^2 yr)","Mass Accumulation Rate Uncertainty u(r(i)) g/(cm^2 yr)","Dry Bulk Density Section ri g/cm^3","Dry Bulk Density Section Uncertainty u(ri) g/cm^3","Dry Bulk Density Layer r(i) g/cm^3","Dry Bulk Density Layer Uncertainty u(r(i)) g/cm^3","Sediment Accumulation Rate s(i) cm/yr","Sediment Accumulation Rate Uncertainty U(s(i)) cm/yr")))
CF_table[,1:2]=mass_table[,1:2]
CF_table[1:nrow(concentrations_table),3]=concentrations_table[,10]
CF_table[1:nrow(concentrations_table),4]=concentrations_table[,11]
CF_table[,12]=mass_table[,4]
CF_table[,13]=mass_table[,7]
#inventory below and uncertainty
for (i in (nrow(CF_table)):1){
  if (i %% 2 ==1){
    if (i==nrow(CF_table)){
      CF_table[nrow(CF_table),5]=missing_inventory
      CF_table[nrow(CF_table),6]=CF_table[nrow(CF_table),5]*sqrt((decay_const_uncer/decay_const)^2+(MAR_final_err/MAR_final)^2+(concentrations_table[nrow(concentrations_table),7]/concentrations_table[nrow(concentrations_table),6])^2)
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
    CF_table[i,11]=CF_table[i,10]*sqrt(((decay_const_uncer/decay_const)^2)+((CF_table[i,6]/CF_table[i,5])^2)+((concentrations_table[i,9]/concentrations_table[i,8])^2))
  }
}
# DBD layer loop
for (i in 1:(nrow(CF_table)-2)){
  if (i %% 2 ==1){
    if (i==1){
      DBD_layer_table=array(data=NA, dim=c(10,3), dimnames=list(NULL, c("zi","mi","DBD g/cm^3")))
      DBD_layer_table[,1]=na.omit(mass_table[1:20,2])
      DBD_layer_table[,2]=na.omit(concentrations_table[1:20,12])
      DBD_layer_table[,3]=na.omit(mass_table[1:20,4])
      DBD_model_list=list()
      for (j in 3:10) {
        # Define the extent of the data to use for the model
        extent_DBD <- j  # This will give 3, 4, ..., 10
        # Extract the x and y values from the array
        x_values_DBD <- DBD_layer_table[1:extent_DBD, 2]  # Column 2 for x values
        y_values_DBD <- DBD_layer_table[1:extent_DBD, 3]  # Column 3 for y values
        # Fit the linear model and store it in the list
        DBD_model_list[[j-2]] <- lm(y_values_DBD ~ x_values_DBD)
      }
      # Step 1: Extract R^2 values
      r_squared_values_DBD <- sapply(DBD_model_list, function(model) summary(model)$adj.r.squared)
      # Step 2: Identify the index of the model with the highest R^2
      best_model_index_DBD <- which.max(r_squared_values_DBD)
      # Step 3: Retrieve the corresponding model
      best_model_DBD <- DBD_model_list[[best_model_index_DBD]]
      CF_table[1,14]=summary(best_model_DBD)$coefficients[1,1]
      #CF_table[i,14]= mean(CF_table[2:6,12], na.rm=TRUE)
    }else{
      CF_table[i,14]=mean(CF_table[(i-1):(i+1),12], na.rm=TRUE)
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
      CF_table[i,15]=summary(best_model_DBD)$coefficients["(Intercept)", "Std. Error"]
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

#save the data
folder_name3 = "CF modeling data"
dir.create(folder_name3)
CF_table_path <- file.path(folder_name3, "CF_table.csv")
write.csv(CF_table, CF_table_path, row.names = FALSE, na = "")
# year plot vs layer depth-z(i)
pdf(file.path(folder_name3, "age_vs_layer_depth.pdf"), width = 8, height = 6)
CF_age_table = array(data = NA, dim = c(length(na.omit(CF_table[, 9])), 3), dimnames = list(NULL, c("x", "x_uncertainty", "y")))
CF_age_table[, 3] = na.omit(CF_table[1:(nrow(CF_table) - 2), 1])
CF_age_table[, 2] = na.omit(CF_table[, 8])
CF_age_table[, 1] = na.omit(CF_table[, 9])
CF_age_table = as.data.frame(CF_age_table)

CF_age_table_plot = ggplot(CF_age_table, aes(x = x, y = y)) +
  geom_point() +  # Add points
  geom_errorbar(aes(xmin = x - x_uncertainty, xmax = x + x_uncertainty), width = 0.2) +  # Add error bars
  labs(title = "Age vs Layer Depth", x = "Year (CE)", y = "Layer Depth (cm)") +  # Add labels
  scale_x_reverse(breaks = seq(min(CF_age_table$x), max(CF_age_table$x), by = 10), 
                  labels = number_format(accuracy = 1)) +  # Set x-axis increments of 10 without decimals
  scale_y_reverse(breaks = seq(min(CF_age_table$y), max(CF_age_table$y), by = interval_thickness*2), labels=scales::number_format(accuracy=0.5)) +  # Set x-axis dynamically
  theme_minimal()  # Use a minimal theme

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

CF_age_md_table_plot = ggplot(CF_age_md_table, aes(x = x, y = y)) +
  geom_point() +  # Add points
  geom_errorbar(aes(xmin = x - x_uncertainty, xmax = x + x_uncertainty), width = 0.2) +  # Add error bars
  labs(title = "Age vs Cumulative Mass Depth", x = "Year (CE)", y = "Cumulative Mass Depth (g/cm^2)") +  # Add labels
  scale_x_reverse(breaks = seq(min(CF_age_md_table$x), max(CF_age_md_table$x), by = 10), 
                  labels = number_format(accuracy = 1)) +  # Set x-axis increments of 10 without decimals
  scale_y_reverse(breaks = seq(min(CF_age_md_table$y), max(CF_age_md_table$y), by = breaks_by_value(mass_depth_table, "Mass Depth mi (g/cm^2)")), labels=scales::number_format(accuracy=accuracy_finder(mass_depth_table[,2]))) +
  theme_minimal()  # Use a minimal theme

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
  geom_point(data = MAR_plot_table, aes(x = x, y=y)) +
  geom_path(data = MAR_plot_table, aes(x=x, y = y), color = "red") +
  geom_errorbar(data = MAR_plot_table,
                aes(y = y, xmin = x - x_uncertainty, xmax = x + x_uncertainty),
                width = 0.2) +
  labs(title = "Mass Accumulation Rate by Year",
       x = "Mass Accumualtion Rate (g/(cm^2 yr))",
       y = "Year (CE)") +
  scale_x_continuous(breaks = seq(min(MAR_plot_table$x), max(MAR_plot_table$x), by = breaks_by_value(MAR_plot_table, "x")), labels=scales::number_format(accuracy=accuracy_finder(MAR_plot_table[,2]))) +
  scale_y_continuous(breaks = seq(min(MAR_plot_table$y), max(MAR_plot_table$y), by = 5),
                     labels = scales::number_format(accuracy = 1)) +
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
  geom_point(data = SAR_plot_table, aes(x = x, y=y)) +
  geom_path(data = SAR_plot_table, aes(x=x, y = y), color = "blue") +
  geom_errorbar(data = SAR_plot_table,
                aes(y = y, xmin = x - x_uncertainty, xmax = x + x_uncertainty),
                width = 0.2) +
  labs(title = "Sediment Accumulation Rate by Year",
       x = "Sediment Accumulation Rate (cm/yr)",
       y = "Year (CE)") +
  scale_x_continuous(breaks = seq(min(SAR_plot_table$x), max(SAR_plot_table$x), by = breaks_by_value(SAR_plot_table,"x")), labels=scales::number_format(accuracy=accuracy_finder(SAR_plot_table[,2]))) +
  scale_y_continuous(breaks = seq(min(SAR_plot_table$y), max(SAR_plot_table$y), by = 5),
                     labels = scales::number_format(accuracy = 1)) +
  theme_minimal()

print(SAR_plot)
dev.off()

model_summaries_DBD <- data.frame(Model = character(),
                                  Coefficient = character(),
                                  Estimate = numeric(),
                                  Std.Error = numeric(),
                                  t.value = numeric(),
                                  Pr.t = numeric(),
                                  Adj.R.squared = numeric(),
                                  stringsAsFactors = FALSE)

# Extract information from each model and store it in the data frame
for (model_index in seq_along(DBD_model_list)) {
  model <- DBD_model_list[[model_index]]
  
  # Check if the model is valid
  if (!is.null(model)) {
    summary_model <- summary(model)
    
    # Extract coefficients
    coefficients <- summary_model$coefficients
    for (i in 1:nrow(coefficients)) {
      model_summaries_DBD <- rbind(model_summaries_DBD, 
                                   data.frame(Model = paste("Model", model_index, "Indices 1-",model_index+2),  # Adjusting index for naming
                                              Coefficient = rownames(coefficients)[i],
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
write_xlsx(model_summaries_DBD, path = "CF modeling data/DBD_model_summaries.xlsx")

# run some CF/CFCS comparisons, save plots and data to their own folder
if (surface_active_zone_check==FALSE){
folder_name5 = "CF_CFCS model comparisons"
dir.create(folder_name5)
CF_CFCS_table=array(data=NA, dim=c(nrow(mass_table),6), dimnames=list(NULL,c("Top of Interval z(i) (cm)", "Mid Depth zi (cm)","Mass Depth Section mi gcm^-2","Cumulative Mass Depth Layer m(i) gcm^-2","CF Age (CE)","CFCS Age (CE)")))
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
CF_comp_z_table=na.omit(as.matrix(CF_CFCS_table[, c("Top of Interval z(i) (cm)", "CF Age (CE)")]))
#CFCS table for depth plot
CFCS_comp_z_table=na.omit(as.matrix(CF_CFCS_table[, c("Mid Depth zi (cm)", "CFCS Age (CE)")]))
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
CF_CFCS_comp_plot_z=ggplot(combined_data, aes(x = Age, y = Depth, color = Source)) +
  geom_point() +
  geom_line(aes(group = Source)) +
  labs(x = "Age (CE)", y = "Depth (cm)", title = "CF and CFCS Age Model Comparisons (z)") +
  theme_minimal() +
  scale_color_manual(values = c("CF" = "blue", "CFCS" = "red")) +  # Customize colors as needed
  scale_x_continuous(trans = "reverse", 
                     limits = c(max(combined_data$Age), min(combined_data$Age)), 
                     breaks = seq(min(combined_data$Age), max(combined_data$Age), by = 10), 
                     labels = number_format(accuracy = 1)) +  # Reverse the x-axis with custom breaks and labels
  scale_y_continuous(trans = "reverse", 
                     limits = c(max(combined_data$Depth), min(combined_data$Depth)), 
                     breaks = seq(min(combined_data$Depth), max(combined_data$Depth), by = 2), 
                     labels = number_format(accuracy = 0.5))  # Reverse the y-axis
print(CF_CFCS_comp_plot_z)
dev.off()
#mass depth comparison plot
# CF table for mass depth plot
CF_comp_md_table=as.matrix(array(data=NA, dim=c((nrow(concentrations_table)/2),2)))
CF_comp_md_table[,1]=na.omit(mass_table[,10])
CF_comp_md_table[,2]=na.omit(CF_table[,9])
#CF_comp_md_table=na.omit(as.matrix(concentrations_table[, c("Mass Depth Section mi gcm^-2", "CF Age (CE)")]))
#CFCS table for mass depth plot
#CFCS_comp_md_table=na.omit(as.matrix(mass_table[, c("Cumulative Mass Depth Layer m(i) gcm^-2", "CFCS Age (CE)")]))
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
CF_CFCS_comp_plot_md=ggplot(combined_data_md, aes(x = Age, y = Depth, color = Source)) +
  geom_point() +
  geom_line(aes(group = Source)) +
  labs(x = "Age (CE)", y = "Cumulative Mass Depth (g/cm^2)", title = "CF and CFCS Age Model Comparisons (md)") +
  theme_minimal() +
  scale_color_manual(values = c("CF" = "blue", "CFCS" = "red")) +  # Customize colors as needed
  scale_x_continuous(trans = "reverse", 
                     limits = c(max(combined_data_md$Age), min(combined_data_md$Age)), 
                     breaks = seq(min(combined_data_md$Age), max(combined_data_md$Age), by = 10), 
                     labels = number_format(accuracy = 1)) +  # Reverse the x-axis with custom breaks and labels
  scale_y_continuous(trans = "reverse", 
                     limits = c(max(combined_data_md$Depth), min(combined_data_md$Depth)), 
                     breaks = seq(min(combined_data_md$Depth), max(combined_data_md$Depth), by = 2), 
                     labels = number_format(accuracy = 0.5))  # Reverse the y-axis
print(CF_CFCS_comp_plot_md)
dev.off()
# save the cf/cfcs table to CSV
write.csv(CF_CFCS_table, "CF_CFCS model comparisons/CF_CFCS_table.csv", row.names = FALSE, na = "")
}
#completion dialogue for CF
finish_CF=FALSE
print("atmospheric Pb-210 flux (dpm cm^2/yr)")
atm_flux=as.numeric(CF_table[1,5]*decay_const)
atm_flux_uncer=atm_flux*sqrt((as.numeric(CF_table[1,6])/as.numeric(CF_table[1,5]))^2+(decay_const_uncer/decay_const)^2)
print(atm_flux)
flux_table=array(data=NA, dim=c(1,2), dimnames=list(NULL, c("Atmospheric Pb-210 Flux (dpm cm^2/yr)", "Atmospheric Pb-210 Flux Uncertainty (dpm cm^2/yr)")))
flux_table[1,1]=atm_flux
flux_table[1,2]=atm_flux_uncer
flux_table_path <- file.path(folder_name3, "atm_pb210_flux.csv")
write.csv(flux_table, flux_table_path, row.names = FALSE, na = "")
while (finish_CF==FALSE){
  print("missing inventory %")
  print(as.numeric((missing_inventory/CF_table[1,5])*100))
  print("surface active zone depth %")
  if (SAZ_depth!=0){
    print(as.numeric((temp_surface_table[SAZ_depth,2]/concentrations_table[nrow(concentrations_table),1])*100))
  }else{
    print(0)
  }
  
  if (as.numeric((missing_inventory/CF_table[1,5])*100)>=20){
    print("missing inventory is greater than 20%! CF model assumptions violated, interpret with caution")
  }
  if (SAZ_depth!=0){
    if (as.numeric((temp_surface_table[SAZ_depth,2]/concentrations_table[nrow(concentrations_table),1])*100)>=20){
      print("surface active zone is greater than 20% of the core! CF model assumptions violated, interpret with caution")
  }
  }
  finish_CF=readline(prompt = "CF modeling complete. Continue to CA? (True) ")
}
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
boundary_table_path <- file.path(folder_name3, "core_boundary_info.csv")
write.csv(core_boundary_info, boundary_table_path, row.names = FALSE, na = "")
#start of CA section
#initialize array
CA_table=array(data=NA, dim=c((nrow(concentrations_table)+1),17), dimnames=list(NULL,c("Mid Depth zi cm", "ln(Ci)","Age ti yr","Age Uncertainty u(ti) yr","CIC Age (CE)","ur(l)","ur(C0)","ur(Ci)","t(i) yr","u(t(i)) yr","delta t yr","Uncertainty delta t yr"," Sediment Accumulation Rate s cm/yr","Sediment Accumulation Rate Uncertainty u(s) cm/yr","Mass Accumualtion Rate r g/(cm^2 yr)","Mass Accumualtion Rate Uncertainty u(r) g/(cm^2 yr)","Mass Depth Section mi gcm^-2")))
# load copied data
CA_table[,1]=mass_table[,2]
CA_table[1:nrow(concentrations_table),2]=concentrations_table[,13]
CA_table[,17]=mass_table[,10]
#start the loops
#age ti yr, CIC year, ur(l), ur(C0), ur(Ci)
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
      CA_table[i,10]=0 ##confirm this change was correct
    }else{
      if (i==nrow(CA_table)){
        CA_table[i,9]=CA_table[(i-1),3]+(CA_table[(i-1),3]-CA_table[(i-3),3])/2
        CA_table[i,10]=CA_table[(i-2),10]
      }else{
        CA_table[i,9]=mean(CA_table[(i+1):(i-1),3], na.rm=TRUE)
        CA_table[i,10]=sqrt((CA_table[(i+1),4])^2+(CA_table[(i-1),4])^2)/2
      }
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
# save off the data
folder_name4 = "CA modeling data"
dir.create(folder_name4)
CA_table_path <- file.path(folder_name4, "CA_table.csv")
write.csv(CA_table, CA_table_path, row.names = FALSE, na = "")

#make the plots
#log actvity/mass depth
pdf(file.path(folder_name4, "excess_log_pb210_md_CA.pdf"), width = 8, height = 6)
elogpb210_md_plotting_table_CA=array(data=NA, dim=c(nrow(CA_table)/2,3), dimnames=list(NULL, c("x","y", "y_uncertainty")))
elogpb210_md_plotting_table_CA[,1]=na.omit(mass_table[,10])
elogpb210_md_plotting_table_CA[,2]=na.omit(concentrations_table[,13])
elogpb210_md_plotting_table_CA[,3]=na.omit(concentrations_table[,14])
elogpb210_md_plotting_table_CA=as.data.frame(elogpb210_md_plotting_table_CA)
elogpb210_md_plot_CA=ggplot(elogpb210_md_plotting_table_CA, aes(x = x, y = y)) +
  geom_point() +  # Add points
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Excess Log-scale Pb-210 vs Mass Depth", x = "Mass Depth (g/cm^2)", y ="Excess Log-scale Pb-210 (dmp/g)") +  # Add labels
  scale_y_continuous(breaks = seq(min(elogpb210_md_plotting_table_CA$y), max(elogpb210_md_plotting_table_CA$y), by = breaks_by_value(elogpb210_md_plotting_table_CA, "y")), labels=scales::number_format(accuracy=accuracy_finder(elogpb210_md_plotting_table_CA[,3]))) +  # Set x-axis dynamically
  scale_x_continuous(breaks = seq(min(elogpb210_md_plotting_table_CA$x), max(elogpb210_md_plotting_table_CA$x), by = breaks_by_value(mass_depth_table, "Mass Depth mi (g/cm^2)")), labels=scales::number_format(accuracy=accuracy_finder(mass_depth_table[,2]))) +  # Set y-axis dynamically
  theme_minimal()  # Use a minimal theme
print(elogpb210_md_plot_CA)
dev.off()
#log actvity/depth
pdf(file.path(folder_name4, "excess_log_pb210_CA.pdf"), width = 8, height = 6)
elogpb210_plotting_table_CA=array(data=NA, dim=c(nrow(CA_table)/2,3), dimnames=list(NULL, c("x","y", "y_uncertainty")))
elogpb210_plotting_table_CA[,1]=na.omit(CA_table[,1])
elogpb210_plotting_table_CA[,2]=na.omit(concentrations_table[,13])
elogpb210_plotting_table_CA[,3]=na.omit(concentrations_table[,14])
elogpb210_plotting_table_CA=as.data.frame(elogpb210_plotting_table_CA)
elogpb210_CA_plot=ggplot(elogpb210_plotting_table_CA, aes(x = x, y = y)) +
  geom_point() +  # Add points
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Excess Log-scale Pb-210 vs Depth", x = "Depth (cm)", y ="Excess Log-scale Pb-210 (dmp/g)") +  # Add labels
  scale_y_continuous(breaks = seq(min(elogpb210_plotting_table_CA$y), max(elogpb210_plotting_table_CA$y), by = breaks_by_value(elogpb210_plotting_table_CA, "y")), labels=scales::number_format(accuracy=accuracy_finder(elogpb210_plotting_table_CA[,3]))) +  # Set x-axis dynamically
  scale_x_continuous(breaks = seq(min(elogpb210_plotting_table_CA$x), max(elogpb210_plotting_table_CA$x), by = interval_thickness*2), labels=scales::number_format(accuracy=0.5)) +
  theme_minimal()  # Use a minimal theme
print(elogpb210_CA_plot)
dev.off()

#save off lms used in plots above
#log actvity/mass depth
elogpb210_md_plotting_table_CA_summary <- data.frame(Model = character(),
                                                     Coefficient = character(),
                                                     Estimate = numeric(),
                                                     Std.Error = numeric(),
                                                     t.value = numeric(),
                                                     Pr.t = numeric(),
                                                     Adj.R.squared = numeric(),
                                                     stringsAsFactors = FALSE)

# Check if the model is valid
elogpb210_md_plotting_table_CA_model = lm(y ~ x, data = elogpb210_md_plotting_table_CA)
summary_model <- summary(elogpb210_md_plotting_table_CA_model)

# Extract coefficients
coefficients <- summary_model$coefficients
for (i in 1:nrow(coefficients)) {
  elogpb210_md_plotting_table_CA_summary <- rbind(elogpb210_md_plotting_table_CA_summary, 
                                                  data.frame(Model = "CA MAR Model",
                                                             Coefficient = rownames(coefficients)[i],
                                                             Estimate = coefficients[i, "Estimate"],
                                                             Std.Error = coefficients[i, "Std. Error"],
                                                             t.value = coefficients[i, "t value"],
                                                             Pr.t = coefficients[i, "Pr(>|t|)"],
                                                             Adj.R.squared = summary_model$adj.r.squared,  # Add R-squared value
                                                             stringsAsFactors = FALSE))
}

# Write the summary to an Excel file
write_xlsx(elogpb210_md_plotting_table_CA_summary, path = "CA modeling data/elogpb210_md_CA_summary.xlsx")
#log actvity/zi depth
elogpb210_plotting_table_CA_summary <- data.frame(Model = character(),
                                                  Coefficient = character(),
                                                  Estimate = numeric(),
                                                  Std.Error = numeric(),
                                                  t.value = numeric(),
                                                  Pr.t = numeric(),
                                                  Adj.R.squared = numeric(),
                                                  stringsAsFactors = FALSE)

# Check if the model is valid
elogpb210_plotting_table_CA_model = lm(y ~ x, data = elogpb210_plotting_table_CA)
summary_model <- summary(elogpb210_plotting_table_CA_model)

# Extract coefficients
coefficients <- summary_model$coefficients
for (i in 1:nrow(coefficients)) {
  elogpb210_plotting_table_CA_summary <- rbind(elogpb210_plotting_table_CA_summary, 
                                               data.frame(Model = "CA SAR model",
                                                          Coefficient = rownames(coefficients)[i],
                                                          Estimate = coefficients[i, "Estimate"],
                                                          Std.Error = coefficients[i, "Std. Error"],
                                                          t.value = coefficients[i, "t value"],
                                                          Pr.t = coefficients[i, "Pr(>|t|)"],
                                                          Adj.R.squared = summary_model$adj.r.squared,  # Add R-squared value
                                                          stringsAsFactors = FALSE))
}

# Write the summary to an Excel file
write_xlsx(elogpb210_plotting_table_CA_summary, path = "CA modeling data/elogpb210_zi_CA_summary.xlsx")

## make CA SAR plot
pdf(file.path(folder_name4, "SAR_plot_CA.pdf"), width = 8, height = 6)
SAR_plot_table_CA = array(data = NA, dim = c((length(na.omit(CA_table[, 13]))), 3), dimnames = list(NULL, c("x", "y", "y_uncertainty")))
SAR_plot_table_CA[, 1] = na.omit(CA_table[1:(length(CA_table[, 13])), 5])
SAR_plot_table_CA[, 2] = na.omit(CA_table[1:(length(CA_table[, 13])), 13])
SAR_plot_table_CA[, 3] = na.omit(CA_table[1:(length(CA_table[, 13])), 14])
SAR_plot_table_CA = as.data.frame(SAR_plot_table_CA)

SAR_plot_CA <- ggplot() +
  geom_point(data = SAR_plot_table_CA, aes(x = x, y = y)) +
  geom_path(data = SAR_plot_table_CA, aes(x = x, y = y), color = "blue") +
  geom_errorbar(data = SAR_plot_table_CA,
                aes(x = x, ymin = y - y_uncertainty, ymax = y + y_uncertainty),
                width = 0.2) +
  labs(title = "Sediment Accumulation Rate by Year",
       x = "Year (CE)",
       y = "Sediment Accumulation Rate cm/yr") +
  scale_x_reverse(breaks = seq(min(SAR_plot_table_CA$x), max(SAR_plot_table_CA$x), by = (length(concentrations_table[,1])*interval_thickness)/5), 
                  labels = scales::number_format(accuracy = 1)) +  # Set x-axis dynamically
  scale_y_continuous(breaks = seq(min(SAR_plot_table_CA$y), max(SAR_plot_table_CA$y), by = breaks_by_value(SAR_plot_table_CA, "y")), 
                     labels = scales::number_format(accuracy = accuracy_finder(SAR_plot_table_CA[, 3]))) +
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
  geom_point(data = MAR_plot_table_CA, aes(x = x, y = y)) +
  geom_path(data = MAR_plot_table_CA, aes(x = x, y = y), color = "red") +
  geom_errorbar(data = MAR_plot_table_CA,
                aes(x = x, ymin = y - y_uncertainty, ymax = y + y_uncertainty),
                width = 0.2) +
  labs(title = "Mass Accumulation Rate by Year",
       x = "Year (CE)",
       y = "Mass Accumulation Rate (g/(cm^2 yr))") +
  scale_x_reverse(breaks = seq(min(MAR_plot_table_CA$x), max(MAR_plot_table_CA$x), by = (length(concentrations_table[,1])*interval_thickness)/5), 
                  labels = scales::number_format(accuracy = 1)) +  # Set x-axis dynamically
  scale_y_continuous(breaks = seq(min(MAR_plot_table_CA$y), max(MAR_plot_table_CA$y), by = breaks_by_value(MAR_plot_table_CA, "y")), labels=scales::number_format(accuracy=accuracy_finder(MAR_plot_table_CA[,3]))) +
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

CA_age_table_plot = ggplot(CA_age_table, aes(x = x, y = y)) +
  geom_point() +  # Add points
  geom_errorbar(aes(xmin = x - x_uncertainty, xmax = x + x_uncertainty), width = 0.2) +  # Add error bars
  labs(title = "Age vs Mid Depth", x = "Year (CE)", y = "Mid Depth (cm)") +  # Add labels
  scale_x_reverse(breaks = seq(min(CA_age_table$x), max(CA_age_table$x), by = 10), 
                  labels = number_format(accuracy = 1)) +  # Set x-axis increments of 10 without decimals
  scale_y_reverse(breaks = seq(min(CA_age_table$y), max(CF_age_table$y), by = interval_thickness*2), labels=scales::number_format(accuracy=0.5)) +  # Set x-axis dynamically
  theme_minimal()  # Use a minimal theme

print(CA_age_table_plot)
dev.off()

#age mass-depth plot
pdf(file.path(folder_name4, "age_vs_mass_depth.pdf"), width = 8, height = 6)
CA_age_md_table = array(data = NA, dim = c(length(na.omit(CA_table[, 5])), 3), dimnames = list(NULL, c("y", "x_uncertainty", "x")))
CA_age_md_table[, 1] = na.omit(CA_table[, 17])
CA_age_md_table[, 2] = na.omit(CA_table[, 4])
CA_age_md_table[, 3] = na.omit(CA_table[, 5])
CA_age_md_table = as.data.frame(CA_age_md_table)

CA_age_md_table_plot = ggplot(CA_age_md_table, aes(x = x, y = y)) +
  geom_point() +  # Add points
  geom_errorbar(aes(xmin = x - x_uncertainty, xmax = x + x_uncertainty), width = 0.2) +  # Add error bars
  labs(title = "Age vs Mass Depth", x = "Year (CE)", y = "Mass Depth (g/cm^2)") +  # Add labels
  scale_x_reverse(breaks = seq(min(CA_age_md_table$x), max(CA_age_md_table$x), by = 10), 
                  labels = number_format(accuracy = 1)) +  # Set x-axis increments of 10 without decimals
  scale_y_reverse(breaks = seq(min(CA_age_md_table$y), max(CA_age_md_table$y), by = breaks_by_value(mass_depth_table, "Mass Depth mi (g/cm^2)")), labels=scales::number_format(accuracy=accuracy_finder(mass_depth_table[,2]))) +
  theme_minimal()  # Use a minimal theme

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
print("modeling complete!")

save.image(file = "pb210analysis_environment.RData")
