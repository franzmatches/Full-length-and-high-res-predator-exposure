rm(list=ls())
##################################################################################################
####################### Extracting data from video output Homalozoon treatment #########################################
##################################################################################################
# ---- 1. Packages ----

require(tidyverse) 
require(data.table)

# ---- 2. Data import and handling----

#set the video resolution to re-scale from pixels to mm 
#pixels
camera_resolution = c(5440, 3060)

video_duration<-12

#millilmeters 
field_of_view = c(49, 27.5)

#scales
xscale<-camera_resolution[1]/field_of_view[1]

yscale<-camera_resolution[2]/field_of_view[2]

area_scale<-(camera_resolution[1]*camera_resolution[2])/(field_of_view[1]*field_of_view[2])

#first and last name of the species in the model
first_sp<-"HOMver"
last_sp<-"PARcau"

#create a list with all the file names (each object that ends in ".txt" in your folder)
files <- fs::dir_ls(path = "Data/ComTrack/homalozoon", glob = "*.txt")

#with rbind we will paste together the resulting data frame that the following code (the "function") creates for each ".txt" file (x)

summary_data_final<-rbindlist(lapply(as.list(files), function(x){
  
  #x <- as.list(files)[[2]]
  #x <- "hom_25_2_15.txt"
  data <- read.table(x, header=TRUE, dec=".", fill = TRUE) #load the first data sheet
  #create video duration variable to use later to calculate speed
  tracking_duration<-video_duration*(max(data$Frame)-min(data$Frame))/max(data$Frame)
  #we split the name of the file (x) to obtain the treatment name, the replicate number, and the patch number (make sure that the elements are separated by "_")
  split_name <- unlist(strsplit(x, split=c("_", "."))) 
  
  ##if we want to add date to have a singular dataframe with the time series we label like this:
  species_name<- split_name[[1]] #take the Xth element before the "_"
  treatment_name <- gsub(".txt", "", split_name[4]) #take the Xth element before the "_"and drop the ".txt"
  time_point<-split_name[[2]]
  replicate<-split_name[[3]]
  
  #remove false detection (i.e. objects that were detected in less than 2 frames throughout the video, we discard them)
  filter_threshold_frames <- 3 #the threshold of detection across the video (e.g. 3 frames in this case)
  
  #remove false detection
  sorted_data <- data %>% group_by(ID) %>% filter(n()> filter_threshold_frames) %>% ungroup()
  #if condition to continue the loop 
  if (nrow(sorted_data) <= 2) { #we want more than two rows in the data frame
    
    print(x)
    print("No individuals tracked in this video")
    
    sorted_data<-sorted_data %>% add_row() %>% mutate_all(~replace(., is.na(.), 0)) %>% 
      select(!(all_of(first_sp):all_of(last_sp)))
    
    drop_patterns <- c("Box", "Label", "Angle")  #all the variables we want to discard
    
    sorted_data <- sorted_data %>% select(!contains(drop_patterns)) #create clean data without the drop patterns
    
    data_final<-sorted_data %>% mutate(species_name = species_name,
                                       treatment = treatment_name,
                                       time_point = time_point,
                                       replicate = replicate
    )
    
    return(data_final)
    
  } else {
    
    print(x)
    print("sorted data has individuals")
    
    ##we select only the columns that we are gonna use, discarding the others
    drop_patterns <- c("Box", "Label", "Angle")  #all the variables we want to discard
    
    sorted_data <- sorted_data %>% select(!contains(drop_patterns)) #create clean data without the drop patterns
    
    #09/02/2023 edit, let's try rerun the data extraction without re-assigning new ID if video flickers
    filtered_gap_ID<-sorted_data
    
    # #manipulate data to obtain a column with species name and a ID_probability column (easiest for analysis)
    dataz<-filtered_gap_ID %>% group_by(Frame,ID) %>% #group the data 
      tidyr::pivot_longer(all_of(first_sp):all_of(last_sp), names_to = "Species",values_to = "Prob")  #pivot the data
    
    ##data handling to obtain ABUNDANCE estimates with max of ID probability per frame
    
    data_ab<- dataz %>%  ungroup() %>% #is good practive to ungroup before regrouping for new variables!
      mutate(time_s = ((Frame-min(Frame))/(((max(Frame)-min(Frame))/tracking_duration))),  ##create time column with seconds instead of frames
             CentroidX_mm = CentroidX/xscale, CentroidY_mm = CentroidY/yscale, ###scale  measurement of centroids to mm instead of pixels
             Length_um = (Length/xscale)*1000, Width_um = (Width/yscale)*1000,  ##scale protists size to micro meters 
             Area_sqrd_um = (Area/area_scale)*1000000
      ) %>%  
      group_by(ID, Species) %>% #for each ID and each Species...
      mutate(average_p = mean(Prob), median_p = median(Prob)) %>% #create average p and median p column calculated across all frames of the video
      ungroup() %>% 
      group_by(Frame, ID) %>% #for each frame and ID
      filter(average_p == max(average_p)) %>% #we take only the highest average ID probs, we may also use  max(median)
      #filter(median_p == max(median_p))
      ungroup() %>% 
      group_by(Frame, Species) %>% #for each frame and each species... 
      mutate(abundance_frame = n()) %>% #create a column with the count of rows per Frame per species (i.e. how many ID of those species we have per frame)
      ungroup() %>% group_by(Species) %>% 
      #create a column with the max , mean and median abundance per species throughout the video + the same for body size measures
      mutate(max_abundance = max(abundance_frame), mean_abundance = mean(abundance_frame), median_abundance = median(abundance_frame)) %>% 
      ungroup() %>% 
      group_by(ID) %>%#for each ID tracked we want measures of the body size
      mutate( mean_width_um = mean(Width_um), median_width_um = median(Width_um),
              mean_length_um = mean(Length_um), median_length_um = median(Length_um),
              max_length_um = max(Length_um), max_width_um = max(Width_um),
              mean_area_sqrd_um = mean(Area_sqrd_um), median_area_sqrd_um = median(Area_sqrd_um))
    
    ##calculate speed of the protists using data.table
    ##turn data_ab into a data.table format
    data_ab.dt = as.data.table(data_ab)
    setorderv(data_ab.dt, c("ID","time_s")) # order by ID and time_s (formerly Frame)
    
    #operation to calculate speed
    data_ab.dt[ ##use all the dataframe
      , `:=`( #create new columns, in particular...
        #colum deltaY with the difference between the horizontal position of each ID between two time steps
        deltaX = CentroidX_mm - shift(CentroidX_mm, 1L, type = "lag")
        #colum deltaY with the difference between the vertical position of each ID between two time steps
        , deltaY = CentroidY_mm - shift(CentroidY_mm, 1L, type = "lag")
        #diffence in time between two rows
        , deltaTime = time_s - shift(time_s, 1L, type = "lag")
      )
      , ID #do this for each ID
    ][
      , velocity := sqrt(deltaX^2 + deltaY^2)/deltaTime #calculate velocity
    ]
    
    data_ab_velocity<-data_ab.dt %>%ungroup() %>% group_by(ID) %>% #for each ID
      mutate(median_speed = median(velocity, na.rm =T),mean_speed = mean(velocity,na.rm =T)) #create column with mean and median speed throughout the video
    
    #create a table to summarize max , mean and median abundance per species and appen the treatment_name,replicate and Patch_n elements we created before
    data_final<-data_ab_velocity %>% mutate(species_name = species_name,
                                            treatment = paste(treatment_name, collapse = ""),
                                            time_point = time_point,
                                            replicate = replicate
    )
    
    return(data_final)
    
  }
  
}),fill = T)

# ---- 3. Plotting and final cleaning ----

#Plotting to check species abundances, body size and speed extracted from video ############################
#put dataframe in tibble version (more ggplot friendly)
summary_data_final<-as_tibble(summary_data_final) %>%
  group_by(time_point,
           treatment,
           replicate,
           Species) %>% 
  mutate(time_point = as.numeric(time_point)) %>% #time point to be treated as a number
  filter(Species != 'NONpro',
  ) %>% #remove non protist particles 
  ungroup () 


#Plot
ggplot(data = summary_data_final %>% group_by(time_point, treatment, replicate, Species) %>% 
         summarise(max_abundance = unique(max_abundance)),
       aes(x = time_point, y = max_abundance, col = Species))+
  geom_point()+
  geom_line()+
  facet_wrap(~treatment*replicate)+
  scale_y_continuous(limits = c(-1, 35))+
  theme_bw()


#### removing outliers ####
summary_data_percentiles <- summary_data_final %>%
  #for each temperature
  group_by(treatment) %>%
  #adjust time to start from 0 and increase by hours
  mutate(time_point = time_point - 1) %>% ungroup() %>% group_by(treatment, Species) %>% 
  #remove the most extreme values (below 1% and above 99%)
  #only doing this for the variables I'm interested in, and recalculate the mean values
  filter(velocity > quantile(velocity, c(0.01), na.rm = T) & velocity < quantile(velocity, c(0.99), na.rm = T),
         Width_um > quantile(Width_um, c(0.01), na.rm = T) & Width_um < quantile(Width_um, c(0.99), na.rm = T),
         max_length_um > quantile(max_length_um, c(0.01), na.rm = T) & max_length_um < quantile(max_length_um, c(0.99), na.rm = T),
         Area_sqrd_um > quantile(Area_sqrd_um, c(0.01), na.rm = T) & Area_sqrd_um < quantile(Area_sqrd_um, c(0.99), na.rm = T)) %>% 
  ungroup() %>% group_by(treatment,replicate, time_point, Species, ID) %>% 
  mutate(mean_speed = mean(velocity),median_speed = median(velocity),
         mean_width_um = mean(Width_um), median_width_um = median(Width_um),
         mean_length_um = mean(Length_um), median_length_um = median(Length_um),
         max_length_um = max(Length_um), max_width_um = max(Width_um),
         mean_area_sqrd_um = mean(Area_sqrd_um), median_area_sqrd_um = median(Area_sqrd_um)) %>% 
  group_by(treatment,replicate, time_point, Species, Frame) %>% 
  mutate(abundance_frame = n()) %>% group_by(treatment,replicate, time_point, Species) %>% 
  mutate(max_abundance = max(abundance_frame)) 

#clean data to analyze (we select the variables we interested in)
data_homalozoon_clean<-summary_data_percentiles %>% 
  select(!c("CentroidX","CentroidY", "Width", "Length", "Area", "Prob","CentroidX_mm", "CentroidY_mm","average_p","median_p",
            "abundance_frame", "mean_abundance", "median_abundance", "deltaX", "deltaY","deltaTime"))


#plot to see difference from raw data
ggplot(data = data_homalozoon_clean %>% group_by(time_point, treatment, replicate, Species) %>%
         dplyr::summarise(max_abundance_videos = unique(max_abundance))
       ,aes(x = time_point, y = max_abundance_videos, col = Species))+
  geom_point()+
  geom_line()+
  facet_wrap(~treatment*replicate)+
  scale_y_continuous(limits = c(-1, 35))+
  theme_bw()

# ---- 4. Write out data----

saveRDS(data_homalozoon_clean, file = "Data/homalozoon_data_IDs_clean.RDS")
