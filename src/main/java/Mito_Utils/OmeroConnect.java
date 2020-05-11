package Mito_Utils;
/*
 * To the extent possible under law, the OME developers have waived
 * all copyright and related or neighboring rights to this tutorial code.
 *
 * See the CC0 1.0 Universal license for details:
 *     http://creativecommons.org/publicdomain/zero/1.0/
 */

import static Mito_Utils.JDialogOmeroConnect.serverName;
import static Mito_Utils.JDialogOmeroConnect.serverPort;
import static Mito_Utils.JDialogOmeroConnect.userID;
import static Mito_Utils.JDialogOmeroConnect.userPass;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.Roi;
import ij.gui.Line;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.PolygonRoi;
import ij.gui.ShapeRoi;
import ij.gui.TextRoi;
import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import loci.formats.in.DefaultMetadataOptions;
import loci.formats.in.MetadataLevel;
import mcib3d.geom.Object3D;
import mcib3d.image3d.ImageByte;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageShort;
import ome.formats.OMEROMetadataStoreClient;
import ome.formats.importer.ImportCandidates;
import ome.formats.importer.ImportConfig;
import ome.formats.importer.ImportLibrary;
import ome.formats.importer.OMEROWrapper;
import ome.formats.importer.cli.ErrorHandler;
import ome.formats.importer.cli.LoggingImportMonitor;
import ome.formats.importer.targets.ImportTarget;
import ome.model.units.BigResult;
import omero.LockTimeout;
import omero.ServerError;
import omero.api.IQueryPrx;
import omero.api.RawFileStorePrx;
import omero.api.RawPixelsStorePrx;
import omero.cmd.Response;
import omero.gateway.Gateway;
import omero.gateway.LoginCredentials;
import omero.gateway.SecurityContext;
import omero.gateway.exception.DSAccessException;
import omero.gateway.exception.DSOutOfServiceException;
import omero.gateway.facility.*;
import omero.gateway.model.*;
import omero.gateway.rnd.Plane2D;
import omero.log.Logger;
import omero.log.SimpleLogger;
import omero.model.*;
import omero.model.enums.ChecksumAlgorithmSHA1160;
import omero.model.enums.UnitsLength;
import omero.sys.ParametersI;
import java.awt.Polygon;
import java.io.*;
import java.nio.ByteBuffer;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import loci.plugins.in.ImporterOptions;

// encrypt


/**
 * A simple connection to an OMERO server using the Java gateway
 *
 * @author The OME Team
 */
public class OmeroConnect {

    /**
     * Reference to the gateway.
     */
    public static Gateway gateway;
    public static SecurityContext securityContext;

    // options
    private boolean log = false;
    private Logger logger;
   

    
    public OmeroConnect() {
        Logger simpleLogger = new SimpleLogger();
        gateway = new Gateway(simpleLogger);
        //System.out.println("omero connect : " + gateway);
    }

    public static boolean notInExcludeList(String imageName, ArrayList<String> excludeList) {
        boolean ok = true;
        for (String exclude : excludeList) {
            if (imageName.contains(exclude)) {
                ok = false;
            }
        }

        return ok;
    }

    /**
     * Creates a connection, the gateway will take care of the services
     * life-cycle.
     *
     * @param hostname
     * @param port
     * @param userName
     * @param password
     * @return connect ok
     * @throws java.lang.Exception 
     */

    public static boolean connect(String hostname, int port, String userName, String password) throws Exception {
        System.out.println("Connecting to OMERO server " + hostname + ":" + port + " with user " + userName);
        boolean connectSuccess = false;
        //Encrypt encrypt = new Encrypt();
        LoginCredentials cred = new LoginCredentials(userName, password, hostname, port);
        cred.getServer().setHostname(hostname);
        if (port > 0) {
            cred.getServer().setPort(port);
        }
        cred.getUser().setUsername(userName);
        cred.getUser().setPassword(password);
        ExperimenterData user = gateway.connect(cred);
        securityContext = new SecurityContext(user.getGroupId());
        password = "";
        if(gateway.isConnected()) 
            connectSuccess = true;
        else 
           IJ.log("Pb Connection to server " + hostname); 
        return(connectSuccess);
    }
    
    public void setLog(boolean log) {
        this.log = log;
    }
   

    public static long getUserId(String name) throws ExecutionException, DSAccessException, DSOutOfServiceException {
        // FIND ID OF USER
        AdminFacility adminFacility = gateway.getFacility(AdminFacility.class);
        ExperimenterData experimenterData = adminFacility.lookupExperimenter(securityContext, name);
        if (experimenterData == null) {
            IJ.log("Could not add user " + name);
            return 0;
        }
        return experimenterData.getId();
    }

    public long getGroupId(String name, String group) throws ExecutionException, DSAccessException, DSOutOfServiceException {
        // FIND ID OF GROUP
        AdminFacility adminFacility = gateway.getFacility(AdminFacility.class);
        ExperimenterData experimenterData = adminFacility.lookupExperimenter(securityContext, name);
        List<GroupData> groups = experimenterData.getGroups();
        for (GroupData groupData : groups) {
            if (groupData.getName().equalsIgnoreCase(group)) return groupData.getGroupId();
        }

        return 0;
    }

    public long getGroupId(String group) throws ExecutionException, DSAccessException, DSOutOfServiceException {
        // FIND ID OF GROUP
        ExperimenterData experimenterData = gateway.getLoggedInUser();
        List<GroupData> groups = experimenterData.getGroups();
        for (GroupData groupData : groups) {
            if (groupData.getName().equalsIgnoreCase(group)) return groupData.getGroupId();
        }

        return 0;
    }

    /**
     * Makes sure to disconnect to destroy sessions.
     */
    public static void disconnect() {
        gateway.disconnect();
    }

    /**
     * Loads the projects owned by the user currently logged in.
     */
    public void loadProjects() throws Exception {
        long id = 0;
        BrowseFacility browse = gateway.getFacility(BrowseFacility.class);
        Collection<ProjectData> projects = browse.getProjects(securityContext);
        for (ProjectData projectData : projects) {
            IJ.log("project " + projectData.getName() + " " + projectData.getDescription());
            Set<DatasetData> data = projectData.getDatasets();
            for (DatasetData datasetData : data) {
                IJ.log("dataset " + datasetData.getName() + " " + datasetData.getDescription());
                Set<ImageData> images = datasetData.getImages();
                IJ.log("images " + images);
                if (datasetData.getName().contains("FISH 3D")) id = datasetData.getId();
                if (images != null) {
                    for (ImageData imageData : images) {
                        IJ.log("images " + imageData.getName() + " " + imageData.getDescription());
                    }
                }
            }
        }
        BrowseFacility browse2 = gateway.getFacility(BrowseFacility.class);
        System.out.println("**********************" + id);
        Collection<ImageData> images = browse.getImagesForDatasets(securityContext, Arrays.asList(id));
        if (images != null) {
            for (ImageData image : images) {
                System.out.println("images " + image.getName() + " " + image.getDescription());
                PixelsData pixels = image.getDefaultPixels();
                int sizeZ = pixels.getSizeZ(); // The number of z-sections.
                int sizeT = pixels.getSizeT(); // The number of timepoints.
                int sizeC = pixels.getSizeC(); // The number of channels.
                int sizeX = pixels.getSizeX(); // The number of pixels along the X-axis.
                int sizeY = pixels.getSizeY(); // The number of pixels along the Y-axis.
                System.out.println(sizeX + " " + sizeY + " " + sizeZ + " " + sizeC + " " + sizeT);
                long pixelsId = pixels.getId();
                RawDataFacility rdf = gateway.getFacility(RawDataFacility.class);

                ImageHandler handler = new ImageShort("test", sizeX, sizeY, sizeZ);
                Plane2D plane;
                for (int z = 0; z < sizeZ; z++) {
                    System.out.println("Reading plane " + z);
                    plane = rdf.getPlane(securityContext, pixels, z, 0, 0);
                    for (int x = 0; x < sizeX; x++) {
                        for (int y = 0; y < sizeY; y++)
                            handler.setPixelIncrement(x, y, z, (float) plane.getPixelValue(x, y));
                    }
                }
                handler.show();
                RawPixelsStorePrx store = null;
            }
        }
    }

    public static ImageHandler getImageZ(ImageData image, int t, int c, int zmin, int zmax) throws Exception {
        return getImageXYZ(image, t, c, 1, 1, 0, -1, 0, -1, zmin, zmax);
    }

    public ImageHandler getImageBin(ImageData image, int t, int c, int binXY, int binZ) throws Exception {
        return getImageXYZ(image, t, c, binXY, binZ, 0, -1, 0, -1, 0, -1);
    }

    public ImageHandler getImageBinZ(ImageData image, int t, int c, int zmin, int zmax, int binXY, int binZ) throws Exception {
        return getImageXYZ(image, t, c, binXY, binZ, 0, -1, 0, -1, zmin, zmax);
    }

    public ImageHandler getImage(ImageData image, int t, int c) throws Exception {
        return getImageXYZ(image, t, c, 1, 1, 0, -1, 0, -1, 0, -1);
    }

    public static double[] getResolutionImage(ImageData imageData) throws ExecutionException, BigResult {
        RawDataFacility rdf = gateway.getFacility(RawDataFacility.class);
        PixelsData pixels = imageData.getDefaultPixels();
        Length pixelsXY = pixels.getPixelSizeX(UnitsLength.MICROMETER);
        Length pixelsZ = pixels.getPixelSizeZ(UnitsLength.MICROMETER);
        double resXY = 1;
        if (pixelsXY != null) resXY = pixelsXY.getValue();
        double resZ = 1;
        if (pixelsZ != null) resZ = pixelsZ.getValue();

        return new double[]{resXY, resXY, resZ};
    }

    public boolean setResolutionImageUM(ImageData imageData, double resXY, double resZ) {
        try {
            // get the pixels
            PixelsData pixels = imageData.getDefaultPixels();
            pixels.setPixelSizeX(new LengthI(resXY, UnitsLength.MICROMETER));
            pixels.setPixelSizeY(new LengthI(resXY, UnitsLength.MICROMETER));
            pixels.setPixelSizeZ(new LengthI(resZ, UnitsLength.MICROMETER));
            // update database
            gateway.getUpdateService(securityContext).saveObject(pixels.asIObject());
        } catch (ServerError serverError) {
            serverError.printStackTrace();
            return false;
        } catch (DSOutOfServiceException e) {
            e.printStackTrace();
            return false;
        }

        return true;
    }

    /**
     * Load the image from OMERO
     *
     * @param image The ImageData information
     * @param t     frame, starts at 1
     * @param c     channel, starts at 1
     * @param binXY binning value in X-Y, 1,2, ...
     * @param binZ  binning value in Z, 1,2,...
     * @param XMin  for crop area, coordinate x of left crop area
     * @param XMax  for crop area, coordinate x of right crop area
     * @param YMin  for crop area, coordinate y of top  crop area
     * @param YMax  for crop area, coordinate of bottom  crop area
     * @param ZMin  for crop area, coordinate of first slice
     * @param ZMax  for crop area, coordinate of last slice
     * @return a cropped, binned iamge from the imageData information
     * @throws Exception
     */
    public static ImageHandler getImageXYZ(ImageData image, int t, int c, int binXY, int binZ, int XMin, int XMax, int YMin, int YMax, int ZMin, int ZMax) throws Exception {
        RawDataFacility rdf = gateway.getFacility(RawDataFacility.class);
        PixelsData pixels = image.getDefaultPixels();
        int sizeZ = pixels.getSizeZ(); // The number of z-sections.
        int sizeT = pixels.getSizeT(); // The number of timepoints.
        int sizeC = pixels.getSizeC(); // The number of channels.
        int sizeX = pixels.getSizeX(); // The number of pixels along the X-axis.
        int sizeY = pixels.getSizeY(); // The number of pixels along the Y-axis.
        Length pixelsXY = pixels.getPixelSizeX(UnitsLength.MICROMETER);
        Length pixelsZ = pixels.getPixelSizeZ(UnitsLength.MICROMETER);

        String pixelType = pixels.getPixelType();
        double resXY = 1;
        if (pixelsXY != null) resXY = pixelsXY.getValue();
        double resZ = 1;
        if (pixelsZ != null) resZ = pixelsZ.getValue();
        // check t and c
        t = Math.min(t, sizeT);
        t = Math.max(t, 1);
        c = Math.min(c, sizeC);
        c = Math.max(c, 1);
        if (XMax < 0) XMax = sizeX;
        if (YMax < 0) YMax = sizeY;
        if (ZMax < 0) ZMax = sizeZ;
        double sxfull = XMax - XMin;
        double syfull = YMax - YMin;
        double szfull = ZMax - ZMin;
        ImageHandler handler;
        if (pixelType.equals("uint16"))
            handler = new ImageShort("test", (int) Math.ceil(sxfull / binXY), (int) Math.ceil(syfull / binXY), (int) Math.ceil(szfull / binZ));
        else if (pixelType.equals("uint8"))
            handler = new ImageByte("test", (int) Math.ceil(sxfull / binXY), (int) Math.ceil(syfull / binXY), (int) Math.ceil(szfull / binZ));
        else if (pixelType.equals("float"))
            handler = new ImageFloat("test", (int) Math.ceil(sxfull / binXY), (int) Math.ceil(syfull / binXY), (int) Math.ceil(szfull / binZ));
        else {
            IJ.log("Cannot handle pixel type " + pixelType);
            return null;
        }
        handler.setScale(resXY * binXY, resZ * binZ, "um");
        //IJ.log("Calibration " + resXY + " " + resZ);
        // check plane size
        double maxSizePlane = 256000000; // maximum transfer size for OMERO
        double planeSize = (int) (sxfull * syfull);
        if (pixelType.equals("uint16")) planeSize *= 2;
        if (planeSize > maxSizePlane) {
            // divide plane in two (recursive)
            int y1 = (int) Math.floor((YMin + YMax) / 2.0);
            ImageHandler handler1 = getImageXYZ(image, t, c, binXY, binZ, XMin, XMax, YMin, y1, ZMin, ZMax);
            handler.insert(handler1, 0, 0, 0, false);
            int insert2 = handler1.sizeY;
            handler1 = getImageXYZ(image, t, c, binXY, binZ, XMin, XMax, y1, YMax, ZMin, ZMax);
            handler.insert(handler1, 0, insert2, 0, false);
        } else { // read plane from Omero
            Plane2D plane;
            for (int z = ZMin; z < ZMax; z += binZ) {
                //IJ.log("Reading plane " + z + " " + XMin + " " + YMin + " " + c + " " + t);
                plane = rdf.getTile(securityContext, pixels, z, t - 1, c - 1, XMin, YMin, (int) sxfull, (int) syfull);
                for (int x = XMin; x < XMax; x += binXY) {
                    for (int y = YMin; y < YMax; y += binXY) {
                        int xx = (x - XMin) / binXY;
                        int yy = (y - YMin) / binXY;
                        int zz = (z - ZMin) / binZ;
                        handler.setPixel(xx, yy, zz, (float) plane.getPixelValue(x - XMin, y - YMin));
                    }
                }
            }
        }

        return handler;
    }


    public static DatasetData findDataset(String name, ProjectData project, boolean strict) {
        Set<DatasetData> data = project.getDatasets();
        for (DatasetData datasetData : data) {
            //IJ.log("testing "+datasetData.getName());
            if ((strict) && (datasetData.getName().equals(name))) {
                //IJ.log("equals "+name);
                return datasetData;
            }
            if ((!strict) && (datasetData.getName().contains(name))) {
                //IJ.log("contains "+name);
                return datasetData;
            }
        }

        return null;
    }

    public static ArrayList<DatasetData> findDatasets(ProjectData project) {
        Set<DatasetData> data = project.getDatasets();
        ArrayList<DatasetData> datasetData = new ArrayList<>(data.size());
        for (DatasetData set : data) {
            datasetData.add(set);
        }

        return datasetData;
    }


    private static ImageData findOneImage(DatasetData data, String name, boolean strict) throws Exception {
        BrowseFacility browse = gateway.getFacility(BrowseFacility.class);
        Collection<ImageData> images = browse.getImagesForDatasets(securityContext, Arrays.asList(data.getId()));
        for (ImageData imageData : images) {
            if ((strict) && (imageData.getName().equals(name))) return imageData;
            if ((!strict) && (imageData.getName().contains(name))) return imageData;
        }

        return null;
    }

    public static ArrayList<ImageData> findAllImages(DatasetData data) throws Exception {
        ArrayList<ImageData> imageList = new ArrayList<>();
        BrowseFacility browse = gateway.getFacility(BrowseFacility.class);
        Collection<ImageData> images = browse.getImagesForDatasets(securityContext, Arrays.asList(data.getId()));
        for (ImageData imageData : images) {
            imageList.add(imageData);
        }

        return imageList;
    }


    public ArrayList<ImageData> findAllImagesExclude(DatasetData data, ArrayList<String> excludeList) throws
            Exception {
        ArrayList<ImageData> imageList = findAllImages(data);
        ArrayList<ImageData> list = notInExcludeList(imageList, excludeList);
        return list;
    }

    public static ArrayList<ImageData> findImagesContainsName(DatasetData data, String name, boolean strict) throws
            Exception {
        ArrayList<ImageData> imageList = new ArrayList<>();
        BrowseFacility browse = gateway.getFacility(BrowseFacility.class);
        Collection<ImageData> images = browse.getImagesForDatasets(securityContext, Arrays.asList(data.getId()));
        for (ImageData imageData : images) {
            if (strict) {
                if (imageData.getName().equals(name))
                    imageList.add(imageData);
            } else {
                if (imageData.getName().contains(name))
                    imageList.add(imageData);
            }
        }

        return imageList;
    }

    public ArrayList<ImageData> findImagesContainsNameExclude(DatasetData data, String
            name, ArrayList<String> excludeList, boolean strict) throws Exception {
        ArrayList<ImageData> imageList = findImagesContainsName(data, name, strict);
        return notInExcludeList(imageList, excludeList);
    }

    private ArrayList<ImageData> notInExcludeList(ArrayList<ImageData> images, ArrayList<String> excludeList) {
        ArrayList<ImageData> list = new ArrayList<>();
        for (ImageData imageData : images) {
            String imageName = imageData.getName();
            boolean ok = notInExcludeList(imageName, excludeList);
            if (ok) {
                IJ.log("including " + imageName);
                list.add(imageData);
            } else {
                IJ.log("excluding " + imageName);
            }
        }

        return list;
    }

    public static ProjectData findProject(String name, boolean strict) throws Exception {
        BrowseFacility browse = gateway.getFacility(BrowseFacility.class);
        Collection<ProjectData> projects = browse.getProjects(securityContext);
        for (ProjectData projectData : projects) {
            if ((!strict) && (projectData.getName().contains(name))) return projectData;
            if ((strict) && (projectData.getName().equals(name))) return projectData;
        }

        return null;
    }

    public static ArrayList<ProjectData> findAllProjects() throws Exception {
        ArrayList<ProjectData> list = new ArrayList<>();
        BrowseFacility browse = gateway.getFacility(BrowseFacility.class);
        Collection<ProjectData> projects = browse.getProjects(securityContext);
        for (ProjectData projectData : projects) {
            list.add(projectData);

        }

        return list;
    }
    
    public static ArrayList<ProjectData> findUserProjects(Long expID) throws Exception {
        ArrayList<ProjectData> list = new ArrayList<>();
        BrowseFacility browse = gateway.getFacility(BrowseFacility.class);
        Collection<ProjectData> projects = browse.getProjects(securityContext, expID);
        for (ProjectData proj : projects) {
            System.out.println(proj.getName());
            list.add(proj);
        }
        return list;
    }


    public ImageData findOneImage(String project, String data, String image, boolean strict) throws Exception {
        ProjectData projectData;
        if (log) logger.info(null, "Looking for project : " + project);
        projectData = findProject(project, strict);
        if (projectData == null) return null;
        if (log) logger.info(projectData, "Found project : " + projectData.getName());
        DatasetData datasetData = findDataset(data, projectData, strict);
        if (datasetData == null) return null;
        if (log) logger.info(datasetData, "Found dataset : " + datasetData.getName());
        ImageData imageData = findOneImage(datasetData, image, strict);
        if (imageData == null) return null;
        if (log) logger.info(imageData, "Found image : " + imageData.getName());

        return imageData;
    }


    public void addRoi(ImageData image) throws Exception {
        DataManagerFacility dm = gateway.getFacility(DataManagerFacility.class);
        ROIFacility roifac = gateway.getFacility(ROIFacility.class);

        //To retrieve the image see above.
        ROIData data = new ROIData();
        data.setImage(image.asImage());
        //Create a rectangle.
        RectangleData rectangle = new RectangleData(10, 10, 10, 10);
        rectangle.setZ(0);
        rectangle.setT(0);
        data.addShapeData(rectangle);
        Collection<ROIData> roiDatas = roifac.saveROIs(securityContext, image.getId(), Arrays.asList(data));
    }
    
    
    /**
    Retrieve image rois
     * @param image
     * @return roiList
     * @throws java.util.concurrent.ExecutionException
     * @throws omero.gateway.exception.DSOutOfServiceException
     * @throws omero.gateway.exception.DSAccessException
    */
    
    public static List<Roi> getImageRois (ImageData image) throws ExecutionException, DSOutOfServiceException, DSAccessException {
        List<Roi> roiList = new ArrayList<>();
        ROIFacility roifac = gateway.getFacility(ROIFacility.class);
        //Retrieve the roi linked to an image
        List<ROIResult> rois_result = roifac.loadROIs(securityContext, image.getId());
        Roi roi = null;
        for (ROIResult rois : rois_result) {
            for (ROIData r : rois.getROIs()) {
                for (ShapeData s : r.getShapes()) {
                    long id = s.getId();
                    int time = s.getT();
                    int z = s.getZ();
                    if (s instanceof RectangleData) {
                        RectangleData rectData = (RectangleData) s;
                        roi = new Roi(rectData.getX(), rectData.getY(), rectData.getWidth(), rectData.getHeight());
                    }
                    if (s instanceof EllipseData) {
                        EllipseData ellipse = (EllipseData) s;
                        roi = new OvalRoi(ellipse.getX(), ellipse.getY(), ellipse.getRadiusX(), ellipse.getRadiusY());
                    }
                    if (s instanceof PointData) {
                        PointData pt = (PointData) s;
                        roi = new PointRoi(pt.getX(), pt.getY());
                    }
                    if (s instanceof LineData) {
                        LineData lineData = (LineData) s;
                        roi = new Line(lineData.getX1(), lineData.getY1(), lineData.getX2(), lineData.getY2());
                    }
                    if (s instanceof MaskData) {
                        MaskData maskData = (MaskData) s;
                        roi = new Roi(maskData.getX(), maskData.getY(), maskData.getWidth(), maskData.getHeight());
                        
                    }
                    if (s instanceof PolygonData) {
                        PolygonData polyData = (PolygonData) s;
                        List<Point2D.Double> points = polyData.getPoints();
                        Polygon poly = new Polygon();
                        for (Point2D.Double pt : points) {
                            poly.addPoint((int)Math.round(pt.getX()), (int)Math.round(pt.getY()));
                        }
                        roi = new PolygonRoi(poly, Roi.POLYGON);
                    }
                }
                roiList.add((Roi)roi);
            }
        }
        return(roiList);
    }
    
    
    public List<RectangleData> getRectangleROIS(ImageData image) throws Exception {
        DataManagerFacility dm = gateway.getFacility(DataManagerFacility.class);
        ROIFacility roifac = gateway.getFacility(ROIFacility.class);


        //Retrieve the roi linked to an image
        List<ROIResult> roiresults = roifac.loadROIs(securityContext, image.getId());
        ROIResult r = roiresults.iterator().next();
        if (r == null) return null;
        Collection<ROIData> rois = r.getROIs();
        Iterator<ROIData> j = rois.iterator();
        while (j.hasNext()) {
            ROIData roi = j.next();
            System.out.println(" ROI " + roi.getUuid());
            List<ShapeData> list = roi.getShapes(0, 0);
            if (list == null) {
                System.out.println(" NO ROI");
                return null;
            }
            // Do something
            for (ShapeData shapeData : list) {
                System.out.println(" Shape " + shapeData.toString());
                System.out.println(" Shape " + shapeData.getROICoordinate().toString());
            }
        }

        return null;
    }

    /*
    public static boolean processTiles(ArrayList<TapasProcessing> process) throws Exception {
        // find input process
        TapasProcessing input = null;
        for (TapasProcessing processing : process) {
            if (processing instanceof OmeroInput) {
                input = processing;
                break;
            }
        }
        if (input == null) {
            IJ.log("NO INPUT FOUND");
            return false;
        }
        if (process.indexOf(input) > 0) {
            IJ.log("INPUT SHOULD BE FIRST PROCESS");
            return false;
        }
        // info from input
        int tmin, tmax;
        if (!input.getParameter(OmeroInput.FILE).equals("*")) {
            tmin = Integer.parseInt(input.getParameter(OmeroInput.TILEMIN));
            tmax = Integer.parseInt(input.getParameter(OmeroInput.TILEMAX));
        } else {
            tmin = 1;
            tmax = Integer.MAX_VALUE;
        }
        for (int tile = tmin; tile <= tmax; tile++) {
            ImageHandler img = null;
            for (TapasProcessing proc : process) {
                IJ.log(proc.getName());
                img = proc.execute(img, tile);
            }
        }
        return true;
    }*/

    /*
    public ImageHandler processTiles(String project, String dataset, String name, int t, int c, int tileMin, int tileMax, ArrayList<TapasProcessing> process) throws Exception {
        ImageHandler input = null;
        for (int tile = tileMin; tile <= tileMax; tile++) { // max 151
            String name2 = name.replace("?tile?", "" + tile);
            ImageData imageData = this.findImage(project, dataset, name2);

            IJ.log("Reading " + imageData.getName());
            ImageHandler imageHandler = this.getImage(imageData, t, c);
            input = imageHandler.duplicate();
            for (TapasProcessing proc : process) {
                IJ.log(proc.getName());
                input = proc.execute(input, tile);
            }
        }
        return input;
    }*/

    public static void deleteImage(ImageData imageData) throws Exception {
        DataManagerFacility dm = gateway.getFacility(DataManagerFacility.class);
        Response rsp = dm.delete(securityContext, imageData.asIObject()).loop(10, 500);
    }

    public static boolean addImageToDataset(String projectName, String datasetName, String dir, String fileName) throws
            Exception {
        return addImageToDataset(projectName, datasetName, dir, fileName, true);
    }


    public static boolean addImageToDataset(String projectName, String datasetName, String dir, String fileName,
                                     boolean overwrite) throws Exception {
        ProjectData project = findProject(projectName, true);
        System.out.println("Found project : " + project.getName());
        DatasetData datasetData = findDataset(datasetName, project, true);
        System.out.println("Found dataset : " + datasetData.getName());
        // check if file exists
        if (overwrite) {
            ImageData imageData = findOneImage(datasetData, fileName, true);
            if (imageData != null)
                deleteImage(imageData);
        }
        // create import
        File file = new File(dir + fileName);
        String[] paths = new String[]{dir + fileName};
        final ImportConfig config = new ome.formats.importer.ImportConfig();
        config.email.set("");
        config.sendFiles.set(true);
        config.sendReport.set(false);
        config.contOnError.set(false);
        config.debug.set(false);

        config.hostname.set(serverName);
        config.port.set(serverPort);
        config.username.set(userID);
        config.password.set(userPass);

        // use .target instead ??
        config.target.set("Dataset:" + datasetData.getId());
        //config.targetClass.set("omero.model.Dataset");
        //config.targetId.set(datasetData.getId());
        ImportTarget importTarget = config.getTarget();
        System.out.println("Found target : " + importTarget);

// the imported image will go into 'orphaned images' unless
// you specify a particular existing dataset like this:
// config.targetClass.set("omero.model.Dataset");
// config.targetId.set(1L);

        final OMEROMetadataStoreClient store;
        try {
            store = config.createStore();
            //store.logVersionInfo(config.getIniVersionNumber());
            OMEROWrapper reader = new OMEROWrapper(config);
            final ImportLibrary library = new ImportLibrary(store, reader);

            ErrorHandler handler = new ErrorHandler(config);
            library.addObserver(new LoggingImportMonitor());

            final ImportCandidates candidates = new ImportCandidates(reader, paths, handler);
            reader.setMetadataOptions(new DefaultMetadataOptions(MetadataLevel.ALL));

            Thread thread = new Thread(() -> {
                if (library.importCandidates(config, candidates))
                    System.out.println("Import OK");
                else
                    System.out.println("Import PB");
                file.delete();
                store.logout();
            });
            thread.start();
        } catch (Exception e) {
            e.printStackTrace();
        }

        return true;
    }


    public boolean createDataset(String projectName, String dataName, String description) throws Exception {
        DataManagerFacility dm = gateway.getFacility(DataManagerFacility.class);
        ProjectData project = findProject(projectName, true);

        //Using the model object (recommended)
        DatasetData datasetData = new DatasetData();
        datasetData.setName(dataName);
        datasetData.setDescription(description);

        ProjectDatasetLink link = new ProjectDatasetLinkI();
        link.setChild(datasetData.asDataset());
        link.setParent(new ProjectI(project.getId(), false));
        IObject r = dm.saveAndReturnObject(securityContext, link);

        String[] paths = new String[]{"/home/boudier.t/DATA/Danni/2D/Time lapse from 27.03 5.30pm to 28.03.2.00am_Macrophage tracking.tif", "/home/boudier.t/DATA/Danni/2D/Time lapse Macrophage tracking Registered.tif"};


        ImportConfig config = new ome.formats.importer.ImportConfig();
        config.email.set("");
        config.sendFiles.set(true);
        config.sendReport.set(false);
        config.contOnError.set(false);
        config.debug.set(false);

        config.hostname.set(serverName);
        config.port.set(serverPort);
        config.username.set(userID);
        config.password.set(userPass);

// the imported image will go into 'orphaned images' unless
// you specify a particular existing dataset like this:
// config.targetClass.set("core.model.Dataset");
// config.targetId.set(1L);

        OMEROMetadataStoreClient store;
        try {
            store = config.createStore();
            //store.logVersionInfo(config.getIniVersionNumber());
            OMEROWrapper reader = new OMEROWrapper(config);
            ImportLibrary library = new ImportLibrary(store, reader);

            ErrorHandler handler = new ErrorHandler(config);
            library.addObserver(new LoggingImportMonitor());

            ImportCandidates candidates = new ImportCandidates(reader, paths, handler);
            reader.setMetadataOptions(new DefaultMetadataOptions(MetadataLevel.ALL));
            if (library.importCandidates(config, candidates))
                System.out.println("Import OK");
            else System.out.println("Import PB");

            store.logout();

        } catch (Exception e) {
            e.printStackTrace();
        }

        return true;
    }

    public void addRoiRectangleToImage(ImageData image, String text, int x0, int y0, int width, int height, int Zmin, int Zmax, int T) throws ExecutionException, DSAccessException, DSOutOfServiceException {
        DataManagerFacility dm = gateway.getFacility(DataManagerFacility.class);
        ROIFacility roifac = gateway.getFacility(ROIFacility.class);

        //To retrieve the image see above.
        ROIData data = new ROIData();
        data.setImage(image.asImage());
        //Create a rectangle.
        for (int Z = Zmin; Z <= Zmax; Z++) {
            RectangleData rectangle = new RectangleData(x0, y0, width, height);
            rectangle.setZ(Z);
            rectangle.setT(T);
            rectangle.setText(text);
            data.addShapeData(rectangle);
        }
        // Save ROI
        ROIData roiData = roifac.saveROIs(securityContext, image.getId(), Arrays.asList(data)).iterator().next();
    }

    public void addRoimask(ImageData image, String text, byte[] mask, int x0, int y0, int width, int height, int Zmin, int Zmax, int T) throws ExecutionException, DSAccessException, DSOutOfServiceException {
        DataManagerFacility dm = gateway.getFacility(DataManagerFacility.class);
        ROIFacility roifac = gateway.getFacility(ROIFacility.class);

        //To retrieve the image see above.
        ROIData data = new ROIData();
        data.setImage(image.asImage());
        //Create a rectangle.
        for (int Z = Zmin; Z <= Zmax; Z++) {
            MaskData maskData = new MaskData(x0, y0, width, height, mask);
            maskData.setZ(Z);
            maskData.setT(T);
            maskData.setText(text);
            data.addShapeData(maskData);
        }
        // Save ROI
        ROIData roiData = roifac.saveROIs(securityContext, image.getId(), Arrays.asList(data)).iterator().next();
    }

    public void addRoimask(ImageData image, Object3D object3D, int C, int T) throws
            ExecutionException, DSAccessException, DSOutOfServiceException {
        ImageHandler seg = object3D.createSegImage(0, 0, 0);
        byte[][] mask = new byte[seg.sizeZ][seg.sizeXY];
        int offX = seg.offsetX;
        int offY = seg.offsetY;
        int offZ = seg.offsetZ;

        for (int z = 0; z < seg.sizeZ; z++) {
            byte[] img = mask[z];
            int c = 0;
            for (int y = 0; y < seg.sizeY; y++) {
                for (int x = 0; x < seg.sizeX; x++) {
                    if (seg.getPixel(x, y, z) > 0)
                        img[c] = (byte) 1;
                    else
                        img[c] = (byte) 0;
                    c++;
                }
            }
        }
        ROIData data = new ROIData();
        data.setImage(image.asImage());
        for (int z = 0; z < seg.sizeZ; z++) {
            MaskData maskData = new MaskData(offX, offY, seg.sizeX, seg.sizeY, mask[z]);
            maskData.setZ(z + offZ);
            maskData.setT(T);
            maskData.setC(C);
            maskData.setText(object3D.getName() + "_" + object3D.getValue());
            data.addShapeData(maskData);
        }
        // Save ROI
        DataManagerFacility dm = gateway.getFacility(DataManagerFacility.class);
        ROIFacility roifac = gateway.getFacility(ROIFacility.class);
        ROIData roiData = roifac.saveROIs(securityContext, image.getId(), Arrays.asList(data)).iterator().next();
    }

    public void addRoimask(ImageData image, String text, byte[] mask, int x0, int y0, int width, int height, int Zmin, int Zmax) throws ExecutionException, DSAccessException, DSOutOfServiceException {
        addRoimask(image, text, mask, x0, y0, width, height, Zmin, Zmax, 0);
    }


    public void addRoiRectangleToImage(ImageData image, String text, int x0, int y0, int width, int height, int Zmin, int Zmax) throws ExecutionException, DSAccessException, DSOutOfServiceException {
        addRoiRectangleToImage(image, text, x0, y0, width, height, Zmin, Zmax, 0);
    }


    public void addKeyValuePairToProject(ProjectData project, HashMap<String, String> map) throws ExecutionException, DSAccessException, DSOutOfServiceException {
        List<NamedValue> result = new ArrayList<NamedValue>();
        for (String key : map.keySet()) {
            result.add(new NamedValue(key, map.get(key)));
        }
        MapAnnotationData data = new MapAnnotationData();
        data.setContent(result);
        data.setDescription("Training Example");
        //Use the following namespace if you want the annotation to be editable
        //in the webclient and insight
        data.setNameSpace(MapAnnotationData.NS_CLIENT_CREATED);
        DataManagerFacility fac = gateway.getFacility(DataManagerFacility.class);
        fac.attachAnnotation(securityContext, data, project);
    }

    public void addKeyValuePairToDataset(DatasetData dataset, HashMap<String, String> map) throws
            ExecutionException, DSAccessException, DSOutOfServiceException {
        List<NamedValue> result = new ArrayList<>();
        for (String key : map.keySet()) {
            result.add(new NamedValue(key, map.get(key)));
        }
        MapAnnotationData data = new MapAnnotationData();
        data.setContent(result);
        data.setDescription("Training Example");
        //Use the following namespace if you want the annotation to be editable
        //in the webclient and insight
        data.setNameSpace(MapAnnotationData.NS_CLIENT_CREATED);
        DataManagerFacility fac = gateway.getFacility(DataManagerFacility.class);
        fac.attachAnnotation(securityContext, data, dataset);
    }

    private void deleteMapAnnotationData(MapAnnotationData data) throws ExecutionException, DSAccessException, DSOutOfServiceException, LockTimeout, InterruptedException {
        DataManagerFacility dm = gateway.getFacility(DataManagerFacility.class);
        Response rsp = dm.delete(securityContext, data.asIObject()).loop(10, 500);
    }

    public void addKeyValuePairToImage(ImageData image, HashMap<String, String> map) throws ExecutionException, DSAccessException, DSOutOfServiceException, LockTimeout, InterruptedException {
        // delete annotation already existing
        Map<String, MapAnnotationData> dataMap = getValuePairs(image, null);
        // check if contains name
        for (String line : dataMap.keySet()) {
            Map<String, String> result = extractKeyValuePairs(line);
            for (String keys : result.keySet()) {
                for (String myMap : map.keySet()) {
                    if (keys.equalsIgnoreCase(myMap)) {
                        //delete map annotation
                        deleteMapAnnotationData(dataMap.get(line));
                    }
                }
            }
        }

        List<NamedValue> namedValueList = new ArrayList<NamedValue>();
        for (String key : map.keySet()) {
            NamedValue namedValue = new NamedValue(key, map.get(key));
            namedValueList.add(namedValue);
        }
        MapAnnotationData annotationData = new MapAnnotationData();
        annotationData.setContent(namedValueList);
        annotationData.setDescription("Key value pair annotation");
        //Use the following namespace if you want the annotation to be editable
        //in the webclient and insight
        annotationData.setNameSpace(MapAnnotationData.NS_CLIENT_CREATED);
        DataManagerFacility fac = gateway.getFacility(DataManagerFacility.class);
        fac.attachAnnotation(securityContext, annotationData, image);
    }


    public Map<String, MapAnnotationData> getValuePairs(ImageData image, ArrayList<String> addUsers) throws DSOutOfServiceException, ExecutionException, DSAccessException {
        //IJ.log("Finding pairs for keyword " + image.getName());
        Map<String, MapAnnotationData> map = new HashMap<>();
        MetadataFacility metadataFacility = gateway.getFacility(MetadataFacility.class);
        // connection
        ExperimenterData user = gateway.getLoggedInUser();
        SecurityContext ctx = new SecurityContext(user.getGroupId());
        // types
        List<Class<? extends AnnotationData>> types = new ArrayList<>();
        types.add(MapAnnotationData.class);
        // users
        long userId = gateway.getLoggedInUser().getId();
        List<Long> userIds = new ArrayList<>();
        userIds.add(userId);
        // check additional user
        //IJ.log("Add users " + addUsers.size());
        if ((addUsers != null) && (!addUsers.isEmpty())) {
            for (String users : addUsers) {
                //IJ.log("Adding user to pairs " + users);
                Long id = getUserId(users);
                if (id > 0) userIds.add(id);
            }
        }
        // get annotations
        List<AnnotationData> list = metadataFacility.getAnnotations(ctx, image, types, userIds);
        if ((list == null) || (list.isEmpty())) IJ.log("No annotations found");
        else {
            for (AnnotationData annotationData : list) {
                if (annotationData instanceof MapAnnotationData) {
                    MapAnnotationData mapAnnotationData = (MapAnnotationData) annotationData;
                    String line = mapAnnotationData.getContentAsString();
                    map.put(line, mapAnnotationData);
                }
            }
        }

        return map;
    }

    public static List<FileAnnotationData> getFileAnnotations(ImageData image, ArrayList<String> addUsers) throws DSOutOfServiceException, ExecutionException, DSAccessException {
        //IJ.log("Finding pairs for keyword " + image.getName());
        ArrayList<FileAnnotationData> files = new ArrayList<>();
        MetadataFacility metadataFacility = gateway.getFacility(MetadataFacility.class);
        // connection
        ExperimenterData user = gateway.getLoggedInUser();
        SecurityContext ctx = new SecurityContext(user.getGroupId());
        // types
        List<Class<? extends AnnotationData>> types = new ArrayList<>();
        types.add(FileAnnotationData.class);
        // users
        long userId = gateway.getLoggedInUser().getId();
        List<Long> userIds = new ArrayList<>();
        userIds.add(userId);
        // check additional user
        //IJ.log("Add users " + addUsers.size());
        if ((addUsers != null) && (!addUsers.isEmpty())) {
            for (String users : addUsers) {
                //IJ.log("Adding user to pairs " + users);
                Long id = getUserId(users);
                if (id > 0) userIds.add(id);
            }
        }
        // get annotations
        List<AnnotationData> list = metadataFacility.getAnnotations(ctx, image, types, userIds);
        if ((list == null) || (list.isEmpty())) IJ.log("No annotations found");
        else {
            for (AnnotationData annotationData : list) {
                if (annotationData instanceof FileAnnotationData) {
                    FileAnnotationData fileAnnotationData = (FileAnnotationData) annotationData;
                    files.add(fileAnnotationData);
                }
            }
        }

        return files;
    }

    public static FileAnnotationData getFileAnnotation(ImageData imageData, String name, ArrayList<String> addUsers) throws ExecutionException, DSAccessException, DSOutOfServiceException {
        ArrayList<FileAnnotationData> list = (ArrayList<FileAnnotationData>) getFileAnnotations(imageData, addUsers);
        for (FileAnnotationData data : list) {
            data.getAttachedFile();
            if (data.getFileName().equalsIgnoreCase(name)) return data;
        }

        return null;
    }

    private OriginalFile getOriginalFile(long id) throws Exception {
        ParametersI param = new ParametersI();
        param.map.put("id", omero.rtypes.rlong(id));
        IQueryPrx svc = gateway.getQueryService(securityContext);
        return (OriginalFile) svc.findByQuery("select p from OriginalFile as p where p.id = :id", param);
    }

    public File readAttachment(FileAnnotationData annotation) throws Exception {
        int INC = 262144;
        FileAnnotationData fa;
        RawFileStorePrx store = gateway.getRawFileService(securityContext);
        File file = File.createTempFile(annotation.getFileName(), ".tmp");

        OriginalFile of;
        try (FileOutputStream stream = new FileOutputStream(file)) {
            fa = annotation;
            //The id of the original file
            of = getOriginalFile(fa.getFileID());
            store.setFileId(fa.getFileID());
            int offset = 0;
            long size = of.getSize().getValue();
            try {
                for (offset = 0; (offset + INC) < size; ) {
                    stream.write(store.read(offset, INC));
                    offset += INC;
                }
            } finally {
                stream.write(store.read(offset, (int) (size - offset)));
            }
        } finally {
            store.close();
        }

        return file;
    }


    public void printAllAnnotations(ImageData image, ArrayList<String> addUsers) throws DSOutOfServiceException, ExecutionException, DSAccessException {
        //IJ.log("Finding pairs for keyword " + image.getName());
        Map<String, String> map = new HashMap<>();
        MetadataFacility metadataFacility = gateway.getFacility(MetadataFacility.class);
        // connection
        ExperimenterData user = gateway.getLoggedInUser();
        SecurityContext ctx = new SecurityContext(user.getGroupId());
        // types
        List<Class<? extends AnnotationData>> types = new ArrayList<>();
        types.add(MapAnnotationData.class);
        types.add(FileAnnotationData.class);
        types.add(TagAnnotationData.class);
        types.add(RatingAnnotationData.class);
        // users
        long userId = gateway.getLoggedInUser().getId();
        List<Long> userIds = new ArrayList<>();
        userIds.add(userId);
        // check additional user
        //IJ.log("Add users " + addUsers.size());
        if ((addUsers != null) && (!addUsers.isEmpty())) {
            for (String users : addUsers) {
                //IJ.log("Adding user to pairs " + users);
                Long id = getUserId(users);
                if (id > 0) userIds.add(id);
            }
        }
        // get annotations
        IJ.log("Looking for annotations");
        List<AnnotationData> list = metadataFacility.getAnnotations(ctx, image, types, userIds);
        IJ.log("Looking for annotations " + list);
        if ((list == null) || (list.isEmpty())) IJ.log("No annotations found");
        else {
            for (AnnotationData annotationData : list) {
                IJ.log("FOUND ANNOTATION " + annotationData.getContentAsString() + " " + annotationData.getClass().getName());
                if (annotationData instanceof MapAnnotationData) {
                    MapAnnotationData mapAnnotationData = (MapAnnotationData) annotationData;
                    // separate with ; first
                    String line = mapAnnotationData.getContentAsString();
                    String[] pairs = line.split(";");
                    for (int p = 0; p < pairs.length; p++) {
                        String[] data = pairs[p].split("=");
                        map.put(data[0], data[1]);
                    }
                }
            }
        }
    }


    public String getValuePair(ImageData image, String name, ArrayList<String> addUsers) throws DSOutOfServiceException, ExecutionException, DSAccessException {
        Map<String, MapAnnotationData> map = getValuePairs(image, addUsers);
        // check if contains name
        for (String line : map.keySet()) {
            if (line.contains(name)) {
                Map<String, String> result = extractKeyValuePairs(line);
                for (String keys : result.keySet()) {
                    if (keys.equalsIgnoreCase(name)) return result.get(keys);
                }
            }
        }

        return null;
    }

    private Map<String, String> extractKeyValuePairs(String line) {
        Map<String, String> map = new HashMap<>();
        String[] data;
        // extract individual pairs
        if (line.contains(";")) {
            data = line.split(";");
        } else data = new String[]{line};
        // extract key value
        for (String test : data) {
            String[] kvs = test.split("=");
            map.put(kvs[0], kvs[1]);
        }

        return map;
    }

    public String getValuePair(ImageData image, String name) throws DSOutOfServiceException, ExecutionException, DSAccessException {
        return getValuePair(image, name, null);
    }


    public void addTagAnnotation(ProjectData project, String text) throws java.util.concurrent.ExecutionException, DSAccessException, DSOutOfServiceException, ServerError {
        DataManagerFacility dm = gateway.getFacility(DataManagerFacility.class);
        TagAnnotationData tagData = new TagAnnotationData(text);
        tagData.setTagDescription(text);
        ProjectAnnotationLink link;
        link = new ProjectAnnotationLinkI();
        link.setChild(tagData.asAnnotation());
        link.setParent(project.asProject());
        link = (ProjectAnnotationLink) dm.saveAndReturnObject(securityContext, link);
    }

    public void addTagAnnotation(DatasetData dataset, String text) throws ExecutionException, DSAccessException, DSOutOfServiceException {
        DataManagerFacility dm = gateway.getFacility(DataManagerFacility.class);
        TagAnnotationData tagData = new TagAnnotationData(text);
        tagData.setTagDescription(text);
        DatasetAnnotationLinkI link;
        link = new DatasetAnnotationLinkI();
        link.setChild(tagData.asAnnotation());
        link.setParent(dataset.asDataset());
        link = (DatasetAnnotationLinkI) dm.saveAndReturnObject(securityContext, link);
    }

    public void addTagAnnotation(ImageData image, String text) throws ExecutionException, DSAccessException, DSOutOfServiceException {
        DataManagerFacility dm = gateway.getFacility(DataManagerFacility.class);
        TagAnnotationData tagData = new TagAnnotationData(text);
        tagData.setTagDescription(text);
        ImageAnnotationLink link;
        link = new ImageAnnotationLinkI();
        link.setChild(tagData.asAnnotation());
        link.setParent(image.asImage());
        link = (ImageAnnotationLink) dm.saveAndReturnObject(securityContext, link);
    }

    public void addFileAnnotation(ImageData image, File file, String fileType) throws ExecutionException, DSAccessException, DSOutOfServiceException {
        addFileAnnotation(image, file, fileType, "attached");
    }


    public static void addFileAnnotation(ImageData image, File file, String fileType, String comment) throws ExecutionException, DSAccessException, DSOutOfServiceException {
        DataManagerFacility dm = gateway.getFacility(DataManagerFacility.class);
        // checking if annotation exists
        FileAnnotationData fileAnnotationData = getFileAnnotation(image, file.getName(), null);
        if (fileAnnotationData != null) {
            IJ.log("The attachment " + fileAnnotationData.getFileName() + " already exists, deleting.");
            dm.delete(securityContext, fileAnnotationData.asIObject());
        }
        IJ.log("Attaching " + file.getAbsolutePath() + " to " + image.getName());
        Future<FileAnnotationData> annotationData = dm.attachFile(securityContext, file, fileType, comment, null, image);
    }

    // TODO check and update
    public static void addFileAnnotation(Dataset dataset, File file) throws java.util.concurrent.ExecutionException, java.io.IOException, DSAccessException, DSOutOfServiceException, ServerError {
        FileAnnotation fa = createFileAnnotation(file);

        //now link the image and the annotation
        DataManagerFacility dm = gateway.getFacility(DataManagerFacility.class);
        DatasetAnnotationLink link = new DatasetAnnotationLinkI();
        link.setChild(fa);
        link.setParent(dataset);
        //save the link back to the server.
        link = (DatasetAnnotationLink) dm.saveAndReturnObject(securityContext, link);
    }

    // TODO check and update
    public static void addFileAnnotation(ProjectData project, File file) throws java.util.concurrent.ExecutionException, java.io.IOException, DSAccessException, DSOutOfServiceException, ServerError {
        FileAnnotation fa = createFileAnnotation(file);

        //now link the image and the annotation
        DataManagerFacility dm = gateway.getFacility(DataManagerFacility.class);
        ProjectAnnotationLink link = new ProjectAnnotationLinkI();
        link.setChild(fa);
        link.setParent(project.asProject());
        //save the link back to the server.
        link = (ProjectAnnotationLink) dm.saveAndReturnObject(securityContext, link);
    }


    private static FileAnnotation createFileAnnotation(File file) throws java.util.concurrent.ExecutionException, java.io.IOException, DSAccessException, DSOutOfServiceException, ServerError {
        int INC = 262144;
        DataManagerFacility dm = gateway.getFacility(DataManagerFacility.class);

        //To retrieve the image see above.
        String name = file.getName();
        String absolutePath = file.getAbsolutePath();
        String path = absolutePath.substring(0, absolutePath.length() - name.length());

        //create the original file object.
        OriginalFile originalFile = new OriginalFileI();
        originalFile.setName(omero.rtypes.rstring(name));
        originalFile.setPath(omero.rtypes.rstring(path));
        originalFile.setSize(omero.rtypes.rlong(file.length()));
        final ChecksumAlgorithm checksumAlgorithm = new ChecksumAlgorithmI();
        checksumAlgorithm.setValue(omero.rtypes.rstring(ChecksumAlgorithmSHA1160.value));
        originalFile.setHasher(checksumAlgorithm);
        originalFile.setMimetype(omero.rtypes.rstring("application/octet-stream")); // or "application/octet-stream"
        //Now we save the originalFile object
        originalFile = (OriginalFile) dm.saveAndReturnObject(securityContext, originalFile);

        //Initialize the service to load the raw data
        RawFileStorePrx rawFileStore = gateway.getRawFileService(securityContext);

        long pos = 0;
        int rlen;
        byte[] buf = new byte[INC];
        ByteBuffer bbuf;
        //Open file and read stream
        try (FileInputStream stream = new FileInputStream(file)) {
            rawFileStore.setFileId(originalFile.getId().getValue());
            while ((rlen = stream.read(buf)) > 0) {
                rawFileStore.write(buf, pos, rlen);
                pos += rlen;
                bbuf = ByteBuffer.wrap(buf);
                bbuf.limit(rlen);
            }
            originalFile = rawFileStore.save();
        } finally {
            rawFileStore.close();
        }
        //now we have an original File in DB and raw data uploaded.
        //We now need to link the Original file to the image using
        //the File annotation object. That's the way to do it.
        FileAnnotation fa = new FileAnnotationI();
        fa.setFile(originalFile);
        fa.setDescription(omero.rtypes.rstring("Excel")); // The description set above e.g. PointsModel
        fa.setNs(omero.rtypes.rstring("results")); // The name space you have set to identify the file annotation.

        //save the file annotation.
        fa = (FileAnnotation) dm.saveAndReturnObject(securityContext, fa);

        return fa;
    }


}