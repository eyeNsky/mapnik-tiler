#!/usr/bin/env python
'''
Based on:
https://trac.openstreetmap.org/browser/applications/rendering/mapnik/generate_tiles_multiprocess.py
With snippet from:
http://www.klokan.cz/projects/gdal2tiles/gdal2tiles.py

The MIT License (MIT)

Copyright (c) 2015 Jon Sellars

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

from math import pi,cos,sin,log,exp,atan
from subprocess import call
import sys, os
import multiprocessing
import sys,getopt
from osgeo import gdal
try:
    import mapnik2 as mapnik
except:
    import mapnik

USAGE = '''
Usage: mapniktiler.py 
    --image-in srcfile           Input image -Required
    --tile-dir /path/to/dir      Output directory -Optional
    --p value                    Number of threads (2) -Optional
    --z-min value                Min Zoom (0) -Optional
    --z-max value                Max Zoom (Based on GSD) -Optional
    --bbox 'xmin ymin xmax ymax' Calculated from extent -Optional

bbox needs quotes!!
The only required argument is --image-in
eg:  --image-in src.tif  --p 4 --z-min 0 --z-max 6 --bbox '-170 15 -52.0 74.0'
'''


try :
    args,byaa = getopt.getopt(sys.argv[1:], '', ['image-in=','tile-dir=','z-min=','z-max=','bbox=','p=',])
    #print byaa
    #print args
    args = dict(args)


    imageIn = args.get('--image-in')
    tile_dir = args.get('--tile-dir')
    zMin = args.get('--z-min')
    zMax = args.get('--z-max')
    p = args.get('--p')
    b_box = args.get('--bbox')

except:
    print USAGE
    #sys.exit()

'''
Start def's
'''

def calcImgExt(img):
    dataset = gdal.Open(img)
    # get epsg code
    try:
        epsg = dataset.GetProjectionRef().split(',')[-1].split('"')[1]
    except:
        epsg = 0
        print dataset.GetDescription(),'has no projection'
    geot = dataset.GetGeoTransform()    
    # Get image height width and heigth in pixels
    rastX = dataset.RasterXSize
    rastY = dataset.RasterYSize
    # Get pixel sizes
    pixelX = geot[1]
    pixelY = geot[5]
    # Get ULX,ULY
    ULX = geot[0] + pixelX/2
    ULY = geot[3] + pixelY/2
    # Calculate LRX,LRY
    LRX = ULX+(pixelX * rastX)
    LRY = ULY+(pixelY * rastY)
    dataset = None

    imgBounds = ULX,ULY,LRX,LRY
    #print imgBounds
    return imgBounds

def calcImgRes(img):
    dataset = gdal.Open(img)
    geot = dataset.GetGeoTransform()    
    # Get pixel size x and assume y is the same
    pixelX = geot[1]

    # release dataset
    dataset = None


    return pixelX

def jsonTemplate():
    jt = '''{
    "description": "", 
    "bounds": "%s,%s,%s,%s", 
    "minzoom": "%s", 
    "version": "1.0.0", 
    "template": "", 
    "maxzoom": "%s", 
    "name": "%s"
}'''
    return jt

def xmlTemplate():
    mapTemp = '''<?xml version="1.0" encoding="utf-8"?>
<!DOCTYPE Map[]>
<Map srs="+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0.0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs +over" background-color="transparent" maximum-extent="-20037508.34,-20037508.34,20037508.34,20037508.34" >




<Style name="bluemarble" filter-mode="first" >
  <Rule>
    <RasterSymbolizer opacity="1" />
  </Rule>
</Style>
<Layer name="bluemarble"
  srs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs">
    <StyleName>bluemarble</StyleName>
    <Datasource>
       <Parameter name="file"><![CDATA[%s]]></Parameter>
       <Parameter name="type"><![CDATA[gdal]]></Parameter>
    </Datasource>
  </Layer>
  
</Map>'''
    return mapTemp

def testMap():
    theTest = '''<!DOCTYPE html>
<html>
<head>

<title>TestMap</title>
<!-- Local
<script type="text/javascript" src="../leaflet/leaflet.js"></script>
<link rel="stylesheet" href="../leaflet/leaflet.css"/>
-->
<link rel="stylesheet" href="http://cdn.leafletjs.com/leaflet-0.7.2/leaflet.css" />
<script src="http://cdn.leafletjs.com/leaflet-0.7.2/leaflet.js"></script>

<style>
body { margin:0; padding:0; }
#map { position:absolute; top:0; bottom:0; width:100%%; }
</style>

</head>

<body>
<div id="map"></div>
<script type='text/javascript'>
 
      var map = L.map("map", {
          center: [%s,%s],
          zoom: %s,
          fadeAnimation: false
      });

      var mqOAM = new L.tileLayer("http://{s}.mqcdn.com/tiles/1.0.0/map/{z}/{x}/{y}.jpg",{
          type: 'base', 
          tileSize:256, 
          minZoom: 0, 
          maxZoom: 18,
          attribution:'Tiles Courtesy of <a href="http://www.mapquest.com/" target="_blank">MapQuest</a> <img src="http://developer.mapquest.com/content/osm/mq_logo.png">',
          subdomains: ['otile1','otile2','otile3','otile4']
      }).addTo(map);
      
      var tileset = new L.tileLayer("{z}/{x}/{y}.png",{
          tileSize:256, 
          minZoom: %s, 
          maxZoom: %s
      }).addTo(map);
 
</script>

</body>

</html>'''
    return theTest ##expects lat,lon,zoom,minZoom ,maxzoom

## Start klokkan snippet ##
initialResolution = 156543.033928041

def Resolution(zoom):
    "Resolution (meters/pixel) for given zoom level (measured at Equator)"
    return initialResolution / (2**zoom)    

def ZoomForPixelSize(pixelSize):
    "Maximal scaledown zoom of the pyramid closest to the pixelSize."
    for i in range(30):
      if pixelSize > Resolution(i):
        return i-1 if i!=0 else 0 # We don't want to scale up
## End klokkan snippet ##

'''
end 'o defs
'''
# handle missing zooms
# only works on decimal degree projections....
if zMin:
    zMin = int(zMin)
else:
    zMin = 0
if zMax:
    zMax = int(zMax)
else:
    dd = calcImgRes(imageIn)
    m = dd*100000
    zMax = ZoomForPixelSize(m)

if b_box:
    
    xMin = b_box.split()[0]
    yMin = b_box.split()[1]
    xMax = b_box.split()[2]
    yMax = b_box.split()[3]

else:
    xMin,yMax,xMax,yMin = calcImgExt(imageIn)
    

bbox = (float(xMin),float(yMin), float(xMax), float(yMax))
#bbox = (-180,-90,180,90)    
print bbox

#zMin = 0
#zMax = 6

### Done gathering cli arguments

imgPath = os.path.abspath(imageIn)
basename = os.path.basename(imageIn).replace('.','_')
## these need a way to get all extensions!!!!!!!!! 
fName,theExt = os.path.splitext(imgPath)
xmlPath = imgPath.replace(theExt,'.xml')
if tile_dir == None:
    tile_dir = imgPath.replace(theExt,'')

# open mapnik xml
xOut = open(xmlPath,'w')

xmlText = xmlTemplate()%imgPath
xOut.write(xmlText)
xOut.close()
theMapFile = xmlPath

## Begin what was mostly 
DEG_TO_RAD = pi/180
RAD_TO_DEG = 180/pi

# Default number of rendering threads to spawn, should be roughly equal to number of CPU cores available
NUM_THREADS = 2
if p:
    NUM_THREADS = int(p)


def minmax (a,b,c):
    a = max(a,b)
    a = min(a,c)
    return a

class GoogleProjection:
    def __init__(self,levels=18):
        self.Bc = []
        self.Cc = []
        self.zc = []
        self.Ac = []
        c = 256
        for d in range(0,levels):
            e = c/2;
            self.Bc.append(c/360.0)
            self.Cc.append(c/(2 * pi))
            self.zc.append((e,e))
            self.Ac.append(c)
            c *= 2
                
    def fromLLtoPixel(self,ll,zoom):
         d = self.zc[zoom]
         e = round(d[0] + ll[0] * self.Bc[zoom])
         f = minmax(sin(DEG_TO_RAD * ll[1]),-0.9999,0.9999)
         g = round(d[1] + 0.5*log((1+f)/(1-f))*-self.Cc[zoom])
         return (e,g)
     
    def fromPixelToLL(self,px,zoom):
         e = self.zc[zoom]
         f = (px[0] - e[0])/self.Bc[zoom]
         g = (px[1] - e[1])/-self.Cc[zoom]
         h = RAD_TO_DEG * ( 2 * atan(exp(g)) - 0.5 * pi)
         return (f,h)



class RenderThread:
    def __init__(self, tile_dir, mapfile, q, printLock, maxZoom):
        self.tile_dir = tile_dir
        self.q = q
        self.mapfile = mapfile
        self.maxZoom = maxZoom
        self.printLock = printLock

    def render_tile(self, tile_uri, x, y, z):
        # Calculate pixel positions of bottom-left & top-right
        p0 = (x * 256, (y + 1) * 256)
        p1 = ((x + 1) * 256, y * 256)

        # Convert to LatLong (EPSG:4326)
        l0 = self.tileproj.fromPixelToLL(p0, z);
        l1 = self.tileproj.fromPixelToLL(p1, z);

        # Convert to map projection (e.g. mercator co-ords EPSG:900913)
        c0 = self.prj.forward(mapnik.Coord(l0[0],l0[1]))
        c1 = self.prj.forward(mapnik.Coord(l1[0],l1[1]))

        # Bounding box for the tile
        if hasattr(mapnik,'mapnik_version') and mapnik.mapnik_version() >= 800:
            bbox = mapnik.Box2d(c0.x,c0.y, c1.x,c1.y)
        else:
            bbox = mapnik.Envelope(c0.x,c0.y, c1.x,c1.y)
        render_size = 256
        self.m.resize(render_size, render_size)
        self.m.zoom_to_box(bbox)
        if(self.m.buffer_size < 128):
            self.m.buffer_size = 128

        # Render image with default Agg renderer
        im = mapnik.Image(render_size, render_size)
        mapnik.render(self.m, im)
        im.save(tile_uri, 'png256')


    def loop(self):
        
        self.m = mapnik.Map(256, 256)
        # Load style XML
        mapnik.load_map(self.m, self.mapfile, True)
        # Obtain <Map> projection
        self.prj = mapnik.Projection(self.m.srs)
        # Projects between tile pixel co-ordinates and LatLong (EPSG:4326)
        self.tileproj = GoogleProjection(self.maxZoom+1)
                
        while True:
            #Fetch a tile from the queue and render it
            r = self.q.get()
            if (r == None):
                self.q.task_done()
                break
            else:
                (name, tile_uri, x, y, z) = r

            exists= ""
            if os.path.isfile(tile_uri):
                exists= "exists"
            else:
                self.render_tile(tile_uri, x, y, z)
            bytes=os.stat(tile_uri)[6]
            empty= ''
            if bytes == 103:
                empty = " Empty Tile "
            self.printLock.acquire()
            #print name, ":", z, x, y, exists, empty
            self.printLock.release()
            self.q.task_done()



def render_tiles(bbox, mapfile, tile_dir, minZoom=1,maxZoom=18, name="unknown", num_threads=NUM_THREADS):
    print "render_tiles(",bbox, mapfile, tile_dir, minZoom,maxZoom, name,")"

    # Launch rendering threads
    queue = multiprocessing.JoinableQueue(32)
    printLock = multiprocessing.Lock()
    renderers = {}
    for i in range(num_threads):
        renderer = RenderThread(tile_dir, mapfile, queue, printLock, maxZoom)
        render_thread = multiprocessing.Process(target=renderer.loop)
        render_thread.start()
        #print "Started render thread %s" % render_thread.getName()
        renderers[i] = render_thread

    if not os.path.isdir(tile_dir):
         os.mkdir(tile_dir)

    gprj = GoogleProjection(maxZoom+1) 

    ll0 = (bbox[0],bbox[3])
    ll1 = (bbox[2],bbox[1])

    for z in range(minZoom,maxZoom + 1):
        px0 = gprj.fromLLtoPixel(ll0,z)
        px1 = gprj.fromLLtoPixel(ll1,z)

        # check if we have directories in place
        zoom = "%s" % z
        if not os.path.isdir(tile_dir + zoom):
            os.mkdir(tile_dir + zoom)
        for x in range(int(px0[0]/256.0),int(px1[0]/256.0)+1):
            # Validate x co-ordinate
            if (x < 0) or (x >= 2**z):
                continue
            # check if we have directories in place
            str_x = "%s" % x
            if not os.path.isdir(tile_dir + zoom + '/' + str_x):
                os.mkdir(tile_dir + zoom + '/' + str_x)
            for y in range(int(px0[1]/256.0),int(px1[1]/256.0)+1):
                # Validate x co-ordinate
                if (y < 0) or (y >= 2**z):
                    continue
                str_y = "%s" % y
                tile_uri = tile_dir + zoom + '/' + str_x + '/' + str_y + '.png'
                # Submit tile to be rendered into the queue
                t = (name, tile_uri, x, y, z)
                queue.put(t)

    # Signal render threads to exit by sending empty request to queue
    for i in range(num_threads):
        queue.put(None)
    # wait for pending rendering jobs to complete
    queue.join()
    for i in range(num_threads):
        renderers[i].join()



if __name__ == "__main__":
	
    home = os.environ['HOME']
    try:
        mapfile = os.environ['MAPNIK_MAP_FILE']
    except KeyError:
        mapfile = theMapFile
    try:
        tile_dir = os.environ['MAPNIK_TILE_DIR']
    except KeyError:
        tile_dir = tile_dir

    if not tile_dir.endswith('/'):
        tile_dir = tile_dir + '/'

    #-------------------------------------------------------------------------
    #
    # Change the following for different bounding boxes and zoom levels
    #
    # Start with an overview
    # World


    # Europe+
    #bbox = (-170, 15, -52.0, 74.0)
    render_tiles(bbox, mapfile, tile_dir, zMin, zMax , basename)
metaTxt = jsonTemplate()%(xMin,yMin,xMax, yMax,zMin,zMax,basename)
metaPath = tile_dir+'/metadata.json'
m = open(metaPath,'w')
m.write(metaTxt)
m.close()

testPath = tile_dir+'/leaflet.html'
testTxt = testMap()% ( ((yMin+yMax)/2),((xMin+xMax)/2),((zMin+zMax)/2),zMin ,zMax)
t = open(testPath,'w')
t.write(testTxt)
t.close()
