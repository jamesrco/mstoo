# mstoo_dev.R
# Created 7/17/17 by Jamie Collins, james.r.collins@aya.yale.edu

# Some basic functions and playground for R-based ms2 methods development

#### some required libraries ####

library(xcms)
library(CAMERA)
library(LOBSTAHS)

library(repmis) # for sourcing files from web ("repmis" is awesome)

#### initial (crude?) functions for evaluation of fragmentation spectra ####

# necessary functions that will allow us to extract the correct ms2 spectra, evaluate transitions, etc

get.ms2Peaklist = function (precursor.index,sample_ID,xcmsRaw.list) {
  
  ms2data.start = xcmsRaw.list[[sample_ID]]@msnScanindex[precursor.index]
  ms2data.end = xcmsRaw.list[[sample_ID]]@msnScanindex[precursor.index+1]-1
  
  scandata =
    data.frame(xcmsRaw.list[[sample_ID]]@env$msnMz[ms2data.start:ms2data.end],
               xcmsRaw.list[[sample_ID]]@env$msnIntensity[ms2data.start:ms2data.end])
  
  colnames(scandata) = c("mz","Intensity")
  
  return(scandata)
  
}

get.topN = function(peaklist,N) {
  
  ordered.peaklist = peaklist[order(peaklist$Intensity, decreasing = TRUE),]
  
  topN.peaklist = ordered.peaklist[1:N,]
  
  return(topN.peaklist)
  
}

eval.PIspecies = function(peaklist,species,ppm,ms2.lookupClasses) {
  
  # check to make sure there are no blank values in the peaklist; if so, excise them
  
  peaklist = peaklist[!is.na(peaklist$mz),]
  
  # retrieve product ion m/z
  
  prod.ion = ms2.lookupClasses$mz_value[
    rownames(ms2.lookupClasses)==species]
  
  if (any(abs((prod.ion-peaklist[,1])/prod.ion*1000000)<ppm)) {
    
    # it's a match
    
    return(1)
    
  } else {
    
    return(0)
    
  }
  
}

eval.CNLspecies = function(peaklist,species,sample_ID,ppm,ms2.lookupClasses) {
  
  # check to make sure there are no blank values in the peaklist; if so, excise them
  
  peaklist = peaklist[!is.na(peaklist$mz),]
  
  # calculate theoretical m/z of the ion that would be produced via the neutral loss
  # throwing in a mean() here in case xcms associated more than one peak with the group in this particular sample
  CNL_product_mz = mean(xcms.peakdata_thisgroup_pos[xcms.peakdata_thisgroup_pos$sample==sample_ID,1])-
    ms2.lookupClasses$mz_value[
      rownames(ms2.lookupClasses)==this.IDclass]
  
  # perform comparison
  
  if (any(abs((CNL_product_mz-peaklist[,1])/CNL_product_mz*1000000)<ppm)) {
    
    # it's a match
    
    return(1)
    
  } else {
    
    return(0)
    
  }
  
}

sumfrag = function(x) {
  
  if (all(is.na(x))) {
    
    NA
    
  } else {
    
    sum(x, na.rm = T)
    
  }
  
}

#### an implementation of the above, from work on some cultures ####

# requires several files: (1) annotated xcmsSet object for data in positive ion mode, (2) LOBSet object in positive ion mode, and (3) list containing positive-mode xcmsRaw objects for all samples generated using: 
#             xraw <- xcmsRaw("yourfile.mzXML", includeMSn=TRUE)
# the routine below assumes the objects are stored in the list in the order in which they appear in the xsAnnotate and LOBSet objecty
# *** this last step in particular needs to be automated for all of this to work reproducibly and predictably

#### prep workspace ####

### load necessary data objects from other Github repo where they reside, using repmis ###

# pre-generated LOBSet 

source_data('https://raw.githubusercontent.com/jamesrco/LipidPhotoOxBox/master/data/nice/Orbi_MS_data/LOBSTAHS_processed/UNC_Marchetti_diatom_cultures_pos_withoddFA_LOBSet.RData')
Marchetti_diatoms_LOBset_pos = Marchetti_diatom_cultures_pos_withoddFA # rename so easier to work with
Marchetti_diatoms_LOBset_posPeaks = getLOBpeaklist(Marchetti_diatoms_LOBset_pos) # generate peaklist

# pull in the parent xsAnnotate objects & list containing the xcmsRaw objects

source_data('https://raw.githubusercontent.com/jamesrco/LipidPhotoOxBox/master/data/nice/Orbi_MS_data/xsAnnotate_objects/UNC_Marchetti_diatom_cultures_pos_withoddFA_xsAnnotate.RData')
Marchetti_diatoms_xsA_pos = UNC_Marchetti_diatom_cultures_pos_withoddFA_xsAnnotate # rename so easier to work with

source_data('https://github.com/jamesrco/LipidPhotoOxBox/raw/master/data/raw/Orbi_MS_data/xcmsRaw_objects/UNC_Marchetti_diatom_cultures_pos_xcmsRaw.RData')
Marchetti_xsR = UNC_Marchetti_diatom_cultures_pos_xcmsRaw

### preallocate three matrices for our results ###

# Marchetti_diatom.detected_pos_ion_fragSpec: to keep track of how many valid ms2 fragmentation spectra were detected for the feature in positive ion mode
# Marchetti_diatom.fragdata_results: number of ms2 spectra in which the PI or CNL criteria were validated

Marchetti_diatom.detected_pos_ion_fragments = as.data.frame(matrix(NA,nrow(Marchetti_diatoms_LOBset_posPeaks),ncol=7))
Marchetti_diatom.fragdata_results = as.data.frame(matrix(NA,nrow(Marchetti_diatoms_LOBset_posPeaks),ncol=7))

### some necessary definitions ###

# first, define types and values of MS fragmentation experiments for each IP-DAG class
# (i.e., constant neutral loss (CNL) or product ion (PI))

# create empty data frame
Marchetti_diatom.frag_lookup_classes = as.data.frame(matrix(NA,8,2))
rownames(Marchetti_diatom.frag_lookup_classes) = 
  c("PG","PE","PC","MGDG","SQDG","DGDG","DGCC","DGTS_DGTA")
colnames(Marchetti_diatom.frag_lookup_classes) =
  c("Frag_exp_type","mz_value")

# now, populate our data frame with necessary values
# per Popendorf et al., Lipids (2013) 48:185–195

Marchetti_diatom.frag_lookup_classes[,1] =
  c("CNL","CNL","PI","CNL","CNL","CNL","PI","PI")
Marchetti_diatom.frag_lookup_classes[,2] =
  c(189.040224,141.019094,184.073321,197.089937,261.051837,359.142212,104.106690,236.149249)

# iterate through the IP-DAG IDs by sample, retrieve necessary data from the xsAnnotate and xcmsRaw objects, evaluate, and record the results

for (i in 1:(nrow(Marchetti_diatoms_LOBset_posPeaks))) {
  # iterate through each LOBSTAHS ID
  
  # retrieve LOBSTAHS compound ID, lipid species
  
  this.ID = Marchetti_diatoms_LOBset_posPeaks$compound_name[i]
  this.IDclass = Marchetti_diatoms_LOBset_posPeaks$species[i]
  
  if (this.IDclass %in% rownames(Marchetti_diatom.frag_lookup_classes)) {
    # an escape if the lipid class isn't accounted for in our input parameter table
    
    # retrieve underlying positive-mode xcms group and peak data
    
    xcms.peakIDs_thisgroup_pos =
      Marchetti_diatoms_xsA_pos@xcmsSet@groupidx[[Marchetti_diatoms_LOBset_posPeaks$xcms_peakgroup[i]]]
    
    xcms.peakdata_thisgroup_pos =
      as.data.frame(Marchetti_diatoms_xsA_pos@xcmsSet@peaks[xcms.peakIDs_thisgroup_pos,])
    
    # now, iterate through all instances of this putatively identified compound that are present in the xcms peakgroup and in one of the samples of interest (i.e., not a QC) 
    
    for (j in 1:nrow(xcms.peakdata_thisgroup_pos)) {
      
      # retrieve the sample ID
      
      samp_ID = xcms.peakdata_thisgroup_pos$sample[j]
      
      # before beginning, check if the peak is a QC; if so, skip to end
      
      if (samp_ID %in% c(1:7)) {
        
        ### first, pull out fragmentation spectra relevant to this peak, if they exist ###
        
        # retrieve raw (uncorrected) retention times for this peak
        
        RT_min_pos.raw = Marchetti_diatoms_xsA_pos@xcmsSet@rt$raw[[samp_ID]][
          Marchetti_diatoms_xsA_pos@xcmsSet@rt$corrected[[samp_ID]]==
            xcms.peakdata_thisgroup_pos$rtmin[j]]
        
        RT_max_pos.raw = Marchetti_diatoms_xsA_pos@xcmsSet@rt$raw[[samp_ID]][
          Marchetti_diatoms_xsA_pos@xcmsSet@rt$corrected[[samp_ID]]==
            xcms.peakdata_thisgroup_pos$rtmax[j]]
        
        RT_ctr_pos.raw = Marchetti_diatoms_xsA_pos@xcmsSet@rt$raw[[samp_ID]][
          Marchetti_diatoms_xsA_pos@xcmsSet@rt$corrected[[samp_ID]]==
            xcms.peakdata_thisgroup_pos$rt[j]]
        
        # retrieve data from ms2 scans which have the correct precursor mz (i.e., the mz of this instance of the LOBSTAHS-ID'd compound (peak) we are currently considering) AND were acquired within the RT window (raw min/max) of the peak
        
        # will use search window of 40(!) ppm
        # confirmed this is necessary based on manual inspection of data
        ms2_lkup_window.ppm = 40
        
        # first, gather possibly relevant precursor scans based strictly on mass difference
        
        # get an index to these scans
        possible.precursors_pos.ind =
          which(abs(Marchetti_xsR[[samp_ID]]@msnPrecursorMz-xcms.peakdata_thisgroup_pos$mz[j])/
                  xcms.peakdata_thisgroup_pos$mz[j]*1000000 < ms2_lkup_window.ppm)
        
        # whittle this list to make sure the scans *also* fall within the rt window for the parent peak
        # will use a little bit of a buffer to capture any scans falling just outside of the RT range
        
        valid.precursors_pos.ind =
          possible.precursors_pos.ind[Marchetti_xsR[[samp_ID]]@msnRt[possible.precursors_pos.ind]>
                                        (RT_min_pos.raw-20) &
                                        Marchetti_xsR[[samp_ID]]@msnRt[possible.precursors_pos.ind]<
                                        (RT_max_pos.raw+20)]
        
        # to get the actual QE scan numbers corresponding to this index, run the below:
        # Marchetti_xsR[[j]]@msnAcquisitionNum[valid.precursors_pos.ind]
        
        # record the number of valid ms2 spectra in Marchetti_diatom.detected_pos_ion_fragments
        # add to record if there's information already recorded from another isomer; otherwise, simply overwrite the NA placeholder
        
        if (length(valid.precursors_pos.ind)>0) {
          
          if (is.na(Marchetti_diatom.detected_pos_ion_fragments[i,samp_ID])) {
            
            Marchetti_diatom.detected_pos_ion_fragments[i,samp_ID] = length(valid.precursors_pos.ind)
            
          } else {
            
            # there's already some data from a previous iteration; add to it
            
            Marchetti_diatom.detected_pos_ion_fragments[i,samp_ID] = 
              Marchetti_diatom.detected_pos_ion_fragments[i,samp_ID] + length(valid.precursors_pos.ind)
            
          }
          
        } else {
          
          # there is no ms2 data for this parent, at least not how we went about finding it
          
          if (is.na(Marchetti_diatom.detected_pos_ion_fragments[i,samp_ID])) {
            
            Marchetti_diatom.detected_pos_ion_fragments[i,samp_ID] = 0
            
          } else {
            
            # there's already some data from a previous iteration; add to it
            
            Marchetti_diatom.detected_pos_ion_fragments[i,samp_ID] = 
              Marchetti_diatom.detected_pos_ion_fragments[i,samp_ID] + 0
            
          }
          
        }
        
        # provide ourselves with an escape here in case there are no valid ms2 scans
        
        if (length(valid.precursors_pos.ind)>0) {
          
          # now, pull out the ms2 data for each of these scans
          # msnScanindex is constructed in such a way that we can recreate the peaklist for
          # a given scan using the following syntax
          
          relevant_ms2data = apply(as.matrix(valid.precursors_pos.ind),1,get.ms2Peaklist,samp_ID,Marchetti_xsR)
          
          ### now, can get onto the business of actually examining the spectra for the diagnostic transitions ###
          
          ### scenario 1: class type is diagnosed via presence of product ion ###
          
          if (Marchetti_diatom.frag_lookup_classes$Frag_exp_type[
            rownames(Marchetti_diatom.frag_lookup_classes)==this.IDclass] == "PI") {
            # this class is diagnosed via presence of a product ion
            
            # apply some logic for product ion scenario: assume that for a product ion-based ID to be "good", we must observe a feature with the mz of the diagnostic ion (+/- some mz tolerance) as one of the top N peaks (by intensity) in at least one of the relevant ms2 scans
            
            # so, extract the top N (right now, 20) most intense features in each scan
            
            top_features.PI = lapply(relevant_ms2data,get.topN,20)
            
            # evaluate: do the list(s) of the top N most intense fragments contain the diagnostic ion?
            
            PI.eval_result = lapply(top_features.PI, eval.PIspecies, species = this.IDclass, ppm = 12,
                                    ms2.lookupClasses = Marchetti_diatom.frag_lookup_classes)
            
            # record result
            
            if (is.na(Marchetti_diatom.fragdata_results[i,samp_ID])) {
              
              Marchetti_diatom.fragdata_results[i,samp_ID] = sum(unlist(PI.eval_result))
              
            } else {
              
              # there's already some data from a previous iteration; add to it
              
              Marchetti_diatom.fragdata_results[i,samp_ID] = 
                Marchetti_diatom.fragdata_results[i,samp_ID] + sum(unlist(PI.eval_result))
              
            }
            
          } else if (Marchetti_diatom.frag_lookup_classes$Frag_exp_type[
            rownames(Marchetti_diatom.frag_lookup_classes)==this.IDclass] == "CNL") {
            
            ### scenario 2: class type is diagnosed via constant neutral loss ###
            
            # assume it's a good ID in this case as long as an ion corresponding to the diagnostic CNL is one of the top N (50, for now) peaks (by intensity) in the + mode ms2 spectrum
            
            top_features.CNL = lapply(relevant_ms2data,get.topN,50)
            
            # evaluate & record
            
            CNL.eval_result = lapply(top_features.CNL, eval.CNLspecies, species = this.IDclass, sample_ID = samp_ID, ppm = 12,
                                     ms2.lookupClasses = Marchetti_diatom.frag_lookup_classes)
            
            if (is.na(Marchetti_diatom.fragdata_results[i,samp_ID])) {
              
              Marchetti_diatom.fragdata_results[i,samp_ID] = sum(unlist(CNL.eval_result))
              
            } else {
              
              # there's already some data from a previous iteration; add to it
              
              Marchetti_diatom.fragdata_results[i,samp_ID] = 
                Marchetti_diatom.fragdata_results[i,samp_ID] + sum(unlist(CNL.eval_result))
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
}