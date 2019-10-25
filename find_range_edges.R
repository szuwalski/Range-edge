#==general code for poor man's range edge analysis w/ kernel densities (both weighted and unweighted)

#==required inputs===========================================================
#  data file of survey data with headings "Year", "Lat", "Lon", "Density"
#  an N x 4 matrix of points that define the axes used for edge analysis 
#  4 columns represent the start lat, start lon, end lat, end lon, respectively
#  N represents the number of axes;  this is user defined
#==============================================================================

#==produces these plots:
#   time series of edge at each axis
#   map with all ranges overlaid (weight or non)
#   folder with all years ranges plotted (weighted or non)

#==needed checks
# is there actually an 'edge' along the defined axis given the surveyed stations? or do the data stop there?
# How to interpret when the range intersects with the axis more than once?
# perhaps guidance on setting levels here?


#==load libraries
library(maps)
library(lattice)
library(PBSmapping)
library(mapdata)    #some additional hires data
library(maptools)   #useful tools such as reading shapefiles
library(mapproj)
library(plotrix)
library(dismo)
library(MASS)
library(RColorBrewer)
library(ggtern)
library(retistruct)
library(ggplot2)

find_edges<-function(data, # data file of survey data with headings "Year", "Lat", "Lon", "Density"
                     axes_ends, # N x 4 matrix with start lat, start lon, end lat, end lon, of 
                     weighted=FALSE, # should the kernel density be weighted by densities or be presence/absence
                     sel_level=0.006, # level of contour to use for the range edge
                     plot_ind_maps=FALSE, # plot yearly maps of distributions? useful for finding good levels for contours
                     input_levels=c(0.0025,0.005,0.01), # controls the levels for individual contour plots, also good for identifying 'sel_level'
                     input_h=2,# bandwidth for kernel density
                     fig_id="name_me") # names for the figures if you have different runs in this same repo
{

 #==initial calcs of kernel density
  years<-na.omit(unique(data$Year))
  unique_stations<-unique(data[,2:3])
  dens2<-list(list())
  for(x in 1:length(years))
  {
    #==subset to a given year and take only complete cases with no NAs
    temp<-data[data$Year==years[x],]
    temp<-temp[complete.cases(temp),]
    
    #==calculate kernel densities
    if(weighted==0)
    {
      dens2[[x]] <- kde2d(temp$Lon, temp$Lat, 
                     lims=c(min(temp$Lon)-sd(temp$Lon), max(temp$Lon)+sd(temp$Lon), min(temp$Lat)-sd(temp$Lat), max(temp$Lat)+sd(temp$Lat) ) ,
                     h=c(bandwidth.nrd(temp$Lon)/input_h, bandwidth.nrd(temp$Lon)/  input_h) ) 
    }
    if(weighted==1)
    {
      dens2[[x]] <- kde2d.weighted(temp$Lon, temp$Lat, w=temp$Density,
                              lims=c(min(temp$Lon)-sd(temp$Lon), max(temp$Lon)+sd(temp$Lon),min(temp$Lat)-sd(temp$Lat), max(temp$Lat)+sd(temp$Lat)),
                              h=c(bandwidth.nrd(temp$Lon)/input_h, bandwidth.nrd(temp$Lon)/  input_h)  )
    }
    
    #==make input levels comparable among years (this could probably be done with 'input_h', actually)
    dens2[[x]]$z<-dens2[[x]]$z/sum(dens2[[x]]$z)
  }  
  
  #==plot map with all densities of level input as 'sel_level'
  #==create directory
  dir.create(file.path("plots"), showWarnings = FALSE)
  png(paste("plots/",fig_id,"MMB_densities.png",sep=""),height=8,width=8,res=1200,units='in')
   cols<-rep(0,length(years))
   cols[1:length(years)]<-colorRampPalette(c("Red","white","Dark Blue"))(length(years))
 
   state.map <- map('worldHires', 
                   xlim=c(min(data$Lon,na.rm=T), max(data$Lon,na.rm=T)), 
                   ylim=c(min(data$Lat,na.rm=T),max(data$Lat,na.rm=T)),
                   plot = TRUE, fill = TRUE,col='grey85')
   points(x=unique_stations[,2],y=unique_stations[,1],pch=16,cex=.4)

  for(x in 1:length(years))
    f<-contour(dens2[[x]],  level=sel_level, add=TRUE,col=cols[x],drawlabels=FALSE) 

   for(x in 1:nrow(axes_ends))
    lines(x=c(axes_ends[x,1],axes_ends[x,3]),y=c(axes_ends[x,2],axes_ends[x,4]),lwd=2,lty=2)
      
   legend("bottomleft",bty='n',col=c(cols[1],cols[length(years)]),pch=15,legend=c(years[1],years[length(years)]))
  dev.off()
  
  
  #==plot individual maps, if indicated
  if(plot_ind_maps)
  {
  #== create if yearly directory does not exist
  dir.create(file.path("plots/year"), showWarnings = FALSE)
    
  cols<-rep(0,length(years))
  cols[1:length(years)]<-colorRampPalette(c("Red","white","Dark Blue"))(length(years))
  
  for(x in 1:length(years))
  {
    temp<-data[data$Year==years[x],]
    temp<-temp[complete.cases(temp),]
    
    #==make maps
    if(weighted==0)
      png(paste("plots/year/MMB_dens_",fig_id,"_",years[x],".png",sep=""),height=8,width=8,res=300,units='in')
    if(weighted==1)
      png(paste("plots/year/MMB_dens_wt",fig_id,"_",years[x],".png",sep=""),height=8,width=8,res=300,units='in')
    state.map <- map('worldHires', 
                     xlim=c(min(data$Lon,na.rm=T), max(data$Lon,na.rm=T)), 
                     ylim=c(min(data$Lat,na.rm=T),max(data$Lat,na.rm=T)),
                     plot = TRUE, fill = TRUE,col='grey85')
    
    points(x=unique_stations[,2],y=unique_stations[,1],pch=16,cex=.4)
    points(x=temp$Lon,y=temp$Lat,col=2,pch=16,cex=.7)
    mtext(side=3,x)
    legend("bottomleft",bty='n',col=c(1,2),pch=c(16),legend=c("Surveyed","Present"))
    
    #==plot contours
    f<-contour(dens2[[x]],  level=input_levels, add=TRUE) 
    dev.off()
  }

  }
  
  #=================================================================
  # calculate range edges given input axes
  #=================================================================
  store_edges<-NULL
  for(y in 1:length(years))
  {
   select_contour<-contourLines(dens2[[y]],levels=sel_level)
   int_temp<-NULL  
   
   #==find the intersection of the range edge and input 'edge axes', store in store_edges
   for(x in 1:nrow(axes_ends))
   {
     int_cnt<-1
    for(z in 1:(length(select_contour[[1]]$x)-1))
    {
     chk_edg<-line.line.intersection(P1=axes_ends[x,1:2],
                                     P2=axes_ends[x,3:4],
                                     P3=c(select_contour[[1]]$x[z],select_contour[[1]]$y[z]),
                                     P4=c(select_contour[[1]]$x[z+1],select_contour[[1]]$y[z+1]),
                                     interior.only=TRUE)
    if(!is.na(chk_edg[1])) 
    {
      int_temp<-rbind(int_temp,c(years[y],x,int_cnt,chk_edg))
      int_cnt<-int_cnt+1
    }
    } 
   }
   store_edges<-rbind(store_edges,int_temp)
   }


  colnames(store_edges)<-c("Year","Axis","Intersection","Lon","Lat")
  indat<-as.data.frame(store_edges)

    
return(indat)
  
}


#====================================================
# simple example
#====================================================
data<-read.csv("snow_dat.csv")
axes_ends<-matrix(c(-159,52.6,-178,62,
                    -178,52.6,-158,63),ncol=4,nrow=2,byrow=T)


edges_out<-find_edges(data,
           axes_ends, 
           weighted=FALSE,
           sel_level=0.0055,
           plot_ind_maps=FALSE,
           input_levels=c(0.0025,0.005,0.01),
           input_h=2,
           fig_id="two_axes")

#===================================================================
# plot the time series of lon of range edge over time for each axis
# this is messy, but you get the idea
#====================================================================

p <- ggplot(edges_out,aes(x=Year,y=Lon)) +
  geom_line() +
  facet_wrap(~Axis+Intersection )
png(paste("plots/example_Lon_series_twoaxes.png",sep=""),height=8,width=8,res=300,units='in')

p +  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.background=element_blank())
dev.off()

#=================================================
# goofy axis example
#=================================================
axes_ends<-matrix(c(-159,53,-169,57.5,
                    -167,56.6,-158,61,
                    -170,58,-164,63),ncol=4,nrow=3,byrow=T)

edges_out<-find_edges(data,
                      axes_ends, 
                      weighted=FALSE,
                      sel_level=0.0055,
                      plot_ind_maps=FALSE,
                      input_levels=c(0.0025,0.005,0.01),
                      input_h=2,
                      fig_id="Short_axes")

png(paste("plots/example_Lon_series_shortaxes.png",sep=""),height=8,width=8,res=300,units='in')

p <- ggplot(edges_out,aes(x=Year,y=Lon)) +
  geom_line() +
  facet_wrap(~Axis+Intersection )

p +  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.background=element_blank())
dev.off()
