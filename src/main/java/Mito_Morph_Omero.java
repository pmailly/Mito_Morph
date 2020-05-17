/*
 * Analyze mitochondrial morphological network
 * 
 * Author Philippe Mailly
 */

import static Mito_Utils.Mito_Processing.computeMitoParameters;
import static Mito_Utils.Mito_Processing.find_nucleus2;
import ij.*;
import ij.io.FileSaver;
import ij.plugin.PlugIn;
import ij.process.AutoThresholder;
import java.io.File;
import java.util.logging.Level;
import java.util.logging.Logger;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageHandler;
import static Mito_Utils.Mito_Processing.flush_close;
import static Mito_Utils.JDialogOmeroConnect.imageData;
import static Mito_Utils.JDialogOmeroConnect.selectedDataset;
import static Mito_Utils.JDialogOmeroConnect.selectedProject;
import static Mito_Utils.Mito_Processing.analyzeSkeleton;
import static Mito_Utils.Mito_Processing.clearOutSide;
import static Mito_Utils.Mito_Processing.getPopFromImage;
import static Mito_Utils.Mito_Processing.median_filter;
import static Mito_Utils.Mito_Processing.threshold;
import Mito_Utils.OmeroConnect;
import static Mito_Utils.OmeroConnect.addFileAnnotation;
import static Mito_Utils.OmeroConnect.gateway;
import static Mito_Utils.OmeroConnect.getImageZ;
import static Mito_Utils.OmeroConnect.getResolutionImage;
import static Mito_Utils.OmeroConnect.securityContext;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.Duplicator;
import ij.plugin.RGBStackMerge;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import ome.model.units.BigResult;
import omero.gateway.exception.DSAccessException;
import omero.gateway.exception.DSOutOfServiceException;
import omero.gateway.facility.MetadataFacility;
import omero.gateway.model.ChannelData;
import omero.gateway.model.ImageData;
import omero.gateway.model.PixelsData;


public class Mito_Morph_Omero implements PlugIn {

    public String tempDir = System.getProperty("java.io.tmpdir")+File.separator;
    public String outFileResults = tempDir+"results.xls";
    private final String imageExt = ".czi";
    public BufferedWriter bufferOutResults;
    public static final Calibration cal = new Calibration();
    // insert results file in lastdataset  image
    private ImageData firstImage;
    
    @Override
    public void run(String arg) {
        try {
            int imageNum = 0;
            ArrayList<String> ch = new ArrayList();
            // Results table
            ResultsTable results = new ResultsTable();
            for (ImageData image : imageData) {
                if (image.getName().endsWith(imageExt)) {
                    try {
                        PixelsData pixels = image.getDefaultPixels();
                        int sizeZ = pixels.getSizeZ();
                        int sizeC = pixels.getSizeC();
                        String rootName = image.getName().replace(imageExt, "");
                        MetadataFacility mf = gateway.getFacility(MetadataFacility.class);
                        String[] channels = new String[sizeC];
                        for(ChannelData chs : mf.getChannelData(securityContext, image.getId())) {
                            channels[chs.getIndex()] = chs.getChannelLabeling();
                        }
                        imageNum++;
                        if (imageNum == 1) {
                            firstImage = image;
                            double[] res = getResolutionImage(image);
                            cal.pixelWidth = res[0];
                            cal.pixelHeight = res[1];
                            cal.pixelDepth = res[2];
                        }
                        
                        // find if roi exist and clearoutside

                        List<Roi> rois = Mito_Utils.OmeroConnect.getImageRois(image);
                        String seriesName = image.getName();
                        /*
                         * Open channels
                         */
                        
                        
                        // DAPI channel ch1
                        System.out.println("Opening Nucleus channel");
                        ImagePlus imgNucOrg = getImageZ(image, 1, 1, 0, sizeZ).getImagePlus();
                        
                        // Open mito channel ch0
                        System.out.println("Opening Mito channel");
                        ImagePlus imgMitoOrg = getImageZ(image, 1, 0, 0, sizeZ).getImagePlus();
                        
                        //If no rois found take the whole image
                        if (rois.isEmpty()) {
                            System.out.println("No roi found, taking the whole image");
                            rois.add(0, new Roi(0, 0, imgNucOrg.getWidth(), imgNucOrg.getHeight()));
                        }
                        
                        // for all rois
                        for( int r = 0; r < rois.size(); r++) {
                            Roi roi = rois.get(r);
                            
                            Objects3DPopulation nucPop = new Objects3DPopulation();
                            new WaitForUserDialog("test").show();
                            ImagePlus imgNuc = new Duplicator().crop(imgNucOrg);
                            
                            nucPop = find_nucleus2(imgNuc, roi);
                            int totalNucPop = nucPop.getNbObjects();
                            System.out.println("Roi " + (r+1)+" Detected nucleus = "+totalNucPop);
                        
                            // Find mitos
                            imgMitoOrg.setRoi(roi);
                            ImagePlus imgMito = new Duplicator().crop(imgMitoOrg);
                            median_filter(imgMito, 1.5);
                            IJ.run(imgMito, "Laplacian of Gaussian", "sigma=4 scale_normalised negate stack");
                            threshold(imgMito, AutoThresholder.Method.RenyiEntropy, false, false);
                            clearOutSide(imgMito, roi);

                            Objects3DPopulation mitoPop = getPopFromImage(imgMito, cal);
                            //objectsSizeFilter(minMito, maxMito, mitoPop, imgMitoOrg, false); 
                            System.out.println("Mito pop = "+ mitoPop.getNbObjects());

                            // Find mito network morphology
                            // Skeletonize
                            double[] skeletonParams = analyzeSkeleton(imgMito, tempDir + rootName);
                                                    
                            // Compute global Mito parameters
                            // nb of mito, mean mito volume, skeleton parameters
                            IJ.showStatus("Writing parameters ...");
                            computeMitoParameters(nucPop.getNbObjects(), mitoPop, skeletonParams, rootName+seriesName, results);
                        
                            // Save objects image
                            ImageHandler imhMitoObjects = ImageHandler.wrap(imgMito).createSameDimensions();
                            imhMitoObjects.setCalibration(imgMitoOrg.getCalibration());
                            ImageHandler imhNucObjects = imhMitoObjects.duplicate();
                            mitoPop.draw(imhMitoObjects, 255);
                            nucPop.draw(imhNucObjects, 255);
                            ImagePlus[] imgColors = {imhMitoObjects.getImagePlus(), null, imhNucObjects.getImagePlus()};
                            ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
                            imgObjects.setCalibration(cal);
                            IJ.run(imgObjects, "Enhance Contrast", "saturated=0.35");
                            FileSaver ImgObjectsFile = new FileSaver(imgObjects);
                            String imageFile = rootName + "_" + seriesName + "_Objects.tif";
                            ImgObjectsFile.saveAsTiff(tempDir + imageFile);
                            // add image to dataset
                            OmeroConnect.addImageToDataset(selectedProject, selectedDataset, tempDir, imageFile);
                        
                            flush_close(imgObjects);
                            flush_close(imhMitoObjects.getImagePlus());
                            flush_close(imhNucObjects.getImagePlus());
                            flush_close(imgNuc);
                            flush_close(imgMito);
                        }
                        flush_close(imgNucOrg);
                        flush_close(imgMitoOrg);
                            
                    } catch (ExecutionException ex) {
                        Logger.getLogger(Mito_Morph_Omero.class.getName()).log(Level.SEVERE, null, ex);
                    } catch (DSOutOfServiceException | DSAccessException | BigResult ex) {
                        Logger.getLogger(Mito_Morph_Omero.class.getName()).log(Level.SEVERE, null, ex);
                    } catch (Exception ex) {
                        Logger.getLogger(Mito_Morph_Omero.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
            }
            if (results.getCounter() > 0) {
                // add table to dataset
                File resultsFile = new File(tempDir+"results.xls");
                results.save(tempDir+"results.xls");
                addFileAnnotation(firstImage, resultsFile, "text/csv", "Results");
            }
            OmeroConnect.disconnect();
            IJ.showStatus("Process done");
        } catch (ExecutionException | DSAccessException | DSOutOfServiceException  ex) {
            Logger.getLogger(Mito_Morph_Omero.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}