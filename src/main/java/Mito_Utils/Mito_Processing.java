package Mito_Utils;


import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.GaussianBlur3D;
import ij.plugin.ZProjector;
import ij.plugin.filter.RankFilters;
import ij.process.AutoThresholder;
import ij.process.ImageProcessor;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.ThreadLocalRandom;
import mcib3d.geom.Object3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.image3d.processing.FastFilters3D;
import mcib3d.image3d.regionGrowing.Watershed3D;
import sc.fiji.analyzeSkeleton.AnalyzeSkeleton_;
import sc.fiji.analyzeSkeleton.SkeletonResult;
        
 /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose dots_Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author phm
 */
public class Mito_Processing {
   
    // Max distance to track object
    public static double maxDist = 5;
    public static boolean watershed = false; 
    // min volume in microns^3 for dots
    private static final double minMito = 0.02;
// max volume in microns^3 for dots
    private static final double maxMito = Double.MAX_VALUE;
    
 /**
 * 
 * @param FileResults
 * @param resultsFileName
 * @param header
 * @return 
 */
public static BufferedWriter writeHeaders(String outFileResults, String header) throws IOException {
    FileWriter FileResults = new FileWriter(outFileResults, false);
    BufferedWriter outPutResults = new BufferedWriter(FileResults); 
    outPutResults.write(header);
    outPutResults.flush();
    return(outPutResults);
}    
    
    /**
     * Dialog 
     */
    public static String dialog() {
        String dir = "";
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.addDirectoryField("Choose Directory Containing Image Files : ", "");
        gd.addCheckbox(" WaterShed split", watershed);
        gd.showDialog();
        dir = gd.getNextString();
        watershed = gd.getNextBoolean();
        return(dir);
    }
    
    /**
     * Find mitos
     * @param imgMitos
     * @param roi 
     * @return  
     */
    public static Objects3DPopulation findMitos(ImagePlus imgMitos, Roi roi) {
        median_filter(imgMitos, 1.5);
        IJ.run(imgMitos, "Laplacian of Gaussian", "sigma=4 scale_normalised negate stack");
        threshold(imgMitos, AutoThresholder.Method.RenyiEntropy, false, false);
        IJ.run(imgMitos, "Options...", "iterations=2 count=1 do=Open stack");
        clearOutSide(imgMitos, roi);
        Objects3DPopulation mitoPop = getPopFromImage(imgMitos, imgMitos.getCalibration());
        objectsSizeFilter(minMito, maxMito, mitoPop, imgMitos, false); 
        return(mitoPop);
    }
        
    
    /*Median filter 
     * 
     * @param img
     * @param size
     */ 
    public static void median_filter(ImagePlus img, double size) {
        RankFilters median = new RankFilters();
        for (int s = 1; s <= img.getNSlices(); s++) {
            img.setZ(s);
            median.rank(img.getProcessor(), size, RankFilters.MEDIAN);
            img.updateAndDraw();
        }
    }
    

     /**
     * Nucleus segmentation 2
     * @param imgNuc
     * @return 
     */
    public static Objects3DPopulation find_nucleus2(ImagePlus imgNuc, Roi roi) {
        ImagePlus img = imgNuc.duplicate();
        median_filter(img, 4);
        ImageStack stack = new ImageStack(imgNuc.getWidth(), imgNuc.getHeight());
        for (int i = 1; i <= img.getStackSize(); i++) {
            img.setZ(i);
            img.updateAndDraw();
            IJ.run(img, "Nuclei Outline", "blur=60 blur2=70 threshold_method=Triangle outlier_radius=15 outlier_threshold=1 max_nucleus_size=500 "
                    + "min_nucleus_size=20 erosion=5 expansion_inner=5 expansion=5 results_overlay");
            img.setZ(1);
            img.updateAndDraw();
            ImagePlus mask = new ImagePlus("mask", img.createRoiMask().getBufferedImage());
            ImageProcessor ip =  mask.getProcessor();
            ip.invertLut();
            for (int n = 0; n < 3; n++) 
                ip.erode();
            stack.addSlice(ip);
        }
        ImagePlus imgStack = new ImagePlus("Nucleus", stack);
        imgStack.setCalibration(imgNuc.getCalibration());
        clearOutSide(imgStack, roi);
        Objects3DPopulation cellPop = new Objects3DPopulation();
        if (watershed) {
            ImagePlus imgWater = WatershedSplit(imgStack, 8);
            imgWater.setCalibration(imgNuc.getCalibration());
            flush_close(imgStack);
            imgWater.setCalibration(imgNuc.getCalibration());
            cellPop = new Objects3DPopulation(imgWater);
            cellPop.removeObjectsTouchingBorders(imgWater, false);
            flush_close(imgWater);
        }
        else {
            cellPop = new Objects3DPopulation(imgStack);
            cellPop.removeObjectsTouchingBorders(imgStack, false);
            flush_close(imgStack);
        }
        flush_close(img);
        return(cellPop);
    }
    
    
    /**
     * Size filter objects
     * remove touching border
     * 
     * @param min
     * @param max
     * @param objects
     * @param img
     * @param border
    */
    public static void objectsSizeFilter(double min, double max, Objects3DPopulation objects, ImagePlus img, boolean border) {
        ImageHandler imh = ImageInt.wrap(img).createSameDimensions();
        Object3D obj = null;
        boolean remove = false;
        if (objects.getNbObjects() > 0) {
            for (int n = 0; n < objects.getNbObjects(); n++) {
                remove = false;
                obj = objects.getObject(n);
                double vol = obj.getVolumeUnit();
                // remove if touching border
                if (border) {
                    if (obj.touchBorders(imh, false)) {
                        remove = true;
                    }
                }
                // Size filter remove small
                if ((vol < min) || (vol > max)) {
                    remove = true;
                }
                if (remove) {
                    objects.removeObject(n);
                    n--;
                }
            }
        }
    }
    
    // Threshold images and fill holes
    public static void threshold(ImagePlus img, AutoThresholder.Method thMed, boolean fill, boolean calcul) {
        //  Threshold and binarize
       String cal = "";
       img.setZ(img.getNSlices()/2);
       img.updateAndDraw();
       IJ.setAutoThreshold(img, thMed.toString()+" dark");
       Prefs.blackBackground = false;
       if (calcul)
           cal = " calculate";
        IJ.run(img, "Convert to Mask", "method="+thMed.toString()+cal+" background=Dark");
        if (fill) {
            IJ.run(img,"Fill Holes", "stack");
        }
    }
    
    
    public static Objects3DPopulation getPopFromImage(ImagePlus img, Calibration cal) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        labels.setCalibration(cal);
        Objects3DPopulation pop = new Objects3DPopulation(labels);
        return pop;
    }
    
    
    // write object labels
    public static void labelsObject (Object3D obj, ImagePlus img, int number, int color) {
        Font tagFont = new Font("SansSerif", Font.PLAIN, 32);
        int[] box = obj.getBoundingBox();
        int z = (int)obj.getCenterZ();
        int x = box[0] - 2;
        int y = box[2] - 2;
        img.setSlice(z+1);
        ImageProcessor ip = img.getProcessor();
        ip.setFont(tagFont);
        ip.setColor(color);
        ip.drawString(Integer.toString(number), x, y);
        img.updateAndDraw();    
    }
    
    // tag object number with random color
    public static void tagsObject(ImageHandler imh, Object3D obj) {        
        int col = ThreadLocalRandom.current().nextInt(2, 255 + 1);
        obj.draw(imh, col);  
    }
    
   /**
     * ramdom color nucleus population
     */
    public static ImagePlus randomColorPop (Objects3DPopulation cellsPop,  ImageHandler img, boolean label) {
        //create image objects population
        img.set332RGBLut();
        img.setCalibration(img.getCalibration());
        for (int i = 0; i < cellsPop.getNbObjects(); i++) {
            Object3D obj = cellsPop.getObject(i);
            int col = ThreadLocalRandom.current().nextInt(2, 255 + 1);
            obj.draw(img, col);
            if (label)
               labelsObject(obj, img.getImagePlus(), (i+1), col); 
        } 
        return(img.getImagePlus());
    } 

 
    public static ImagePlus doZProjection(ImagePlus img) {
        ZProjector zproject = new ZProjector();
        zproject.setMethod(ZProjector.MAX_METHOD);
        zproject.setStartSlice(1);
        zproject.setStopSlice(img.getNSlices());
        zproject.setImage(img);
        zproject.doProjection();
       return(zproject.getProjection());
    }
    
// Flush and close images
    public static void flush_close(ImagePlus img) {
        img.flush();
        img.close();
    } 
    
    
    /**
     * Clear out side roi
     * @param img
     * @param roi
     */
    public static void clearOutSide(ImagePlus img, Roi roi) {
        roi.setLocation(0, 0);
        for (int n = 1; n <= img.getNSlices(); n++) {
            ImageProcessor ip = img.getImageStack().getProcessor(n);
            ip.setRoi(roi);
            ip.setBackgroundValue(0);
            ip.setColor(0);
            ip.fillOutside(roi);
        }
        img.deleteRoi();
        img.updateAndDraw();
    }
    
    public static ImagePlus WatershedSplit(ImagePlus binaryMask, float rad) {
        float resXY = 1;
        float resZ = 1;
        float radXY = rad;
        float radZ = rad;
        Calibration cal = binaryMask.getCalibration();
        if (cal != null) {
            resXY = (float) cal.pixelWidth;
            resZ = (float) cal.pixelDepth;
            radZ = radXY * (resXY / resZ);
        }
        ImageInt imgMask = ImageInt.wrap(binaryMask);
        ImageFloat edt = EDT.run(imgMask, 0, resXY, resZ, false, 0);
        ImageHandler edt16 = edt.convertToShort(true);
        ImagePlus edt16Plus = edt16.getImagePlus();
        GaussianBlur3D.blur(edt16Plus, 2.0, 2.0, 2.0);
        edt16 = ImageInt.wrap(edt16Plus);
        edt16.intersectMask(imgMask);
        // seeds
        ImageHandler seedsImg = FastFilters3D.filterImage(edt16, FastFilters3D.MAXLOCAL, radXY, radXY, radZ, 0, false);
        Watershed3D water = new Watershed3D(edt16, seedsImg, 0, 0);
        water.setLabelSeeds(true);
        return(water.getWatershedImage3D().getImagePlus());
    }


    /**
     * nalayze skeleton
     * @param img
     * @param output
     * @return {#branch, branchLenght, #endPoint, #junction}
     */

    public static double[] analyzeSkeleton (ImagePlus img, String outFileName) {
        IJ.run(img, "Skeletonize (2D/3D)", "");
	String imgTitle = img.getTitle();
        AnalyzeSkeleton_ analyzeSkeleton = new AnalyzeSkeleton_();
        AnalyzeSkeleton_.calculateShortestPath = true;
        analyzeSkeleton.setup("",img);
        SkeletonResult skeletonResults = analyzeSkeleton.run(AnalyzeSkeleton_.NONE, false, true, null, true, false);
        ImageStack imgStackLab = analyzeSkeleton.getLabeledSkeletons();
        //  compute parameters for each skeleton
        IJ.showStatus("Computing parameters for each skeleton ...");
        ImagePlus imgLab = new ImagePlus(imgTitle+"_LabelledSkel.tif", imgStackLab);
        ImagePlus imgLabProj = doZProjection(imgLab);
        IJ.run(imgLabProj, "3-3-2 RGB", "");
        imgLabProj.setCalibration(img.getCalibration());
        flush_close(imgLab);
        int[] branchNumbers = skeletonResults.getBranches();
        double[] branchLengths = skeletonResults.getAverageBranchLength();
        int[] nbEndPoints = skeletonResults.getEndPoints();
        int[] junctions = skeletonResults.getJunctions();
        int branches = 0;
        double branchLength = 0;
        int endPoint = 0;
        int junction = 0;
        for (int i = 0; i < skeletonResults.getGraph().length; i++) {
            branches += branchNumbers[i];
            branchLength += branchLengths[i];
            endPoint += nbEndPoints[i];
            junction += junctions[i];
        }
        double[] skeletonParams = {branches, branchLength, endPoint, junction};
        FileSaver imgSave = new FileSaver(imgLabProj);
        imgSave.saveAsTiff(outFileName+"_LabelledSkel.tif");
        flush_close(imgLabProj); 
        return(skeletonParams);
    }
    
    /**
    * Local compute mito parameters add to results file
    * @param nuc nucleus number 
    * @param mitoPop mito population
    * @param mitoParams branch number, branch lenght, end points, junctions
     * @param imgName
    * @param results buffer
    * @throws java.io.IOException
    **/
    public static void computeMitoParameters(int nuc, Objects3DPopulation mitoPop, double[] mitoParams, int roi, 
            String imgName, BufferedWriter results) throws IOException {
        IJ.showStatus("Computing mitochondria parameters ....");
        // mito volume
        double mitoVol = 0;
        for (int i = 0; i < mitoPop.getNbObjects(); i++) {
            mitoVol += mitoPop.getObject(i).getVolumeUnit();
        }
        results.write(imgName+"\t"+roi+"\t"+nuc+"\t"+mitoPop.getNbObjects()+"\t"+mitoVol+"\t"+mitoParams[0]+"\t"+mitoParams[1]+"\t"+mitoParams[2]+"\t"+
                mitoParams[3]+"\n");
        results.flush();
    }
    
    /**
    * Compute mito parameters add results to data omero table
    * @param nuc nucleus number 
    * @param mitoPop mito population
    * @param mitoParams branch number, branch lenght, end points, junctions
     * @param rootName
     * @param table
    * @throws java.io.IOException
    **/
    public static void computeMitoParameters(int nuc, Objects3DPopulation mitoPop, double[] mitoParams, 
            String rootName, ResultsTable table) throws IOException {
        IJ.showStatus("Computing mitochondria parameters ....");
        // mito volume
        double mitoVol = 0;
        for (int i = 0; i < mitoPop.getNbObjects(); i++) {
            mitoVol += mitoPop.getObject(i).getVolumeUnit();
        }
        // Add results to table 
        int counter = table.getCounter();
        table.setValue("Image Name", counter,rootName);
        table.setValue("Nucleus number", counter,nuc);
        table.setValue("Mito number", counter,mitoPop.getNbObjects());
        table.setValue("Mito Volume", counter,mitoVol);
        table.setValue("Mito branch number", counter,mitoParams[0]);
        table.setValue("Mito branch length", counter,mitoParams[1]);
        table.setValue("Mito end points", counter,mitoParams[2]);
        table.setValue("Mito junction number", counter,mitoParams[3]);
        //table.updateResults();

    }
    
    
}
