/****************************************************************************
Copyright (c) 2010, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package gph;

import java.awt.Color;
import java.nio.IntBuffer;
import java.util.Iterator;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.ogl.*;
import edu.mines.jtk.sgl.*;
import static edu.mines.jtk.ogl.Gl.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A group node for displaying images on oblate spheroids.
 * One such spheroid is the WGS-84 model of the earth's surface.
 * Images are assumed to be uniformly sampled in (longitude,latitude);
 * that is, the simplest equirectangular projection is assumed.
 * Units of longitude and latitude are degrees.
 *
 * @author Dave Hale, Colorado School of Mines.
 * @version 2010.12.10
 */
public class OblateImageGroup extends Group {

  /**
   * Constructs an oblate image group.
   * @param model the oblate model.
   */
  public OblateImageGroup(OblateModel model) {
    _model = model;
  }

  /**
   * Adds a child node for the specified image.
   * @param image the image.
   * @return the node added.
   */
  public Node addImage(SampledImage image) {
    return addImage(image,1.0);
  }

  /**
   * Adds a child node for the specified image and scale factor.
   * The scale factor enables the specified image to appear above 
   * or below other images in this group.
   * @param image the image.
   * @param scale factor by which to scale (x,y,z) coordinates.
   * @return the node added.
   */
  public Node addImage(SampledImage image, double scale) {
    Node node = new ImageNode(_model,image,scale);
    addChild(node);
    return node;
  }
  
  /**
   * Removes the specified node if it is an image node child of this group.
   * @param node the node to remove.
   */
  public void removeImage(Node node) {
    Iterator<Node> ni = getChildren();
    while (ni.hasNext()) {
      Node n = ni.next();
      if (n instanceof ImageNode) {
        ImageNode imageNode = (ImageNode)n;
        if (imageNode==node)
          removeChild(node);
      }
    }
  }

  /**
   * Adds a child lat-lon grid node with default parameters.
   * The default grid sampling intervals are 30 degrees.
   * The default color, width and scale are black, 1.0, and 1.001,
   * respectively.
   * @return the node added.
   */
  public Node addGrid() {
    return addGrid(Color.BLACK,1.0,1.001);
  }

  /**
   * Adds a child lat-lon grid node with default grid sampling.
   * The default lat-lon grid sampling intervals are 30 degrees.
   * @param color color of grid lines.
   * @param width width of grid lines.
   * @param scale factor by which to scale (x,y,z) coordinates.
   * @return the node added.
   */
  public Node addGrid(Color color, double width, double scale) {
    Sampling glon = new Sampling(13,30.0,0.0);
    Sampling glat = new Sampling(7,30.0,-90.0);
    return addGrid(glon,glat,color,width,scale);
  }

  /**
   * Adds a child lat-lon grid node.
   * @param slon sampling of grid lines for longitude.
   * @param slat sampling of grid lines for latitude.
   * @param color color of grid lines.
   * @param width width of grid lines.
   * @param scale factor by which to scale (x,y,z) coordinates.
   * @return the node added.
   */
  public Node addGrid(
    Sampling slon, Sampling slat, Color color, double width, double scale) 
  {
    Node node = new GridNode(_model,slon,slat,scale);
    ColorState cs = new ColorState();
    cs.setColor(color);
    LineState ls = new LineState();
    ls.setWidth((float)width);
    StateSet ss = new StateSet();
    ss.add(cs);
    ss.add(ls);
    node.setStates(ss);
    addChild(node);
    return node;
  }
  
  /**
   * Removes the specified node if it is a grid node child of this group.
   * @param node the node to remove.
   */
  public void removeGrid(Node node) {
    Iterator<Node> ni = getChildren();
    while (ni.hasNext()) {
      Node n = ni.next();
      if (n instanceof GridNode) {
        GridNode gridNode = (GridNode)n;
        if (gridNode==node)
          removeChild(node);
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  // private

  private OblateModel _model;
  private Node _node1;
  private Node _node2;

  private static class ImageNode extends Node {

    public ImageNode(
      OblateModel model, SampledImage image, double scale) 
    {
      _model = model;
      _image = image;
      _scale = scale;
      _bs = new BoundingSphere(0.0,0.0,0.0,_model.getRadius());
    }

    protected BoundingSphere computeBoundingSphere(boolean finite) {
      return _bs;
    }

    protected void draw(DrawContext dc) {

      // If necessary, make display list and texture.
      if (_displayList==null) {
        Sampling sx = _image.getSamplingX();
        Sampling sy = _image.getSamplingY();
        float dx = (float)sx.getDelta();
        float flon = (float)sx.getFirst()-0.5f*dx;
        float llon = (float)sx.getLast() +0.5f*dx;
        //int nlon = 5; // DEBUGGING
        int nlon = 2+(int)((llon-flon)/3.0f); // start with 3 deg interval
        float dlon = (llon-flon)/(nlon-1); // actual lon sampling interval
        float dy = (float)sy.getDelta();
        float flat = (float)sy.getFirst()-0.5f*dy;
        float llat = (float)sy.getLast() +0.5f*dy;
        if (flat<-90.0f) flat = -90.0f;
        if (llat> 90.0f) llat =  90.0f;
        //int nlat = 5; // DEBUGGING
        int nlat = 2+(int)((llat-flat)/3.0f); // start with 3 deg interval
        float dlat = (llat-flat)/(nlat-1); // actual lat sampling interval
        Sampling slon = new Sampling(nlon,dlon,flon); // lon sampling
        Sampling slat = new Sampling(nlat,dlat,flat); // lat sampling
        int w = _image.getWidth2();
        int h = _image.getHeight2();
        IntBuffer p = _image.getIntBuffer();
        glPixelStorei(GL_UNPACK_ALIGNMENT,1);
        _textureName = new GlTextureName();
        glBindTexture(GL_TEXTURE_2D,_textureName.name());
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP);
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP);
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
        glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,w,h,0,GL_RGBA,GL_UNSIGNED_BYTE,p);
        glBindTexture(GL_TEXTURE_2D,0);
        _displayList = new GlDisplayList();
        glNewList(_displayList.list(),GL_COMPILE);
        for (int ilat=0; ilat<nlat-1; ++ilat) {
          float lat1 = (float)slat.getValue(ilat+1);
          float lat2 = (float)slat.getValue(ilat  );
          float t1 = _image.getTextureT(lat1);
          float t2 = _image.getTextureT(lat2);
          glBegin(GL_TRIANGLE_STRIP);
          for (int ilon=0; ilon<nlon; ++ilon) {
            float lon = (float)slon.getValue(ilon);
            float s = _image.getTextureS(lon);
            double[] xyz = _model.xyz(lat1,lon,0.0);
            float x = (float)xyz[0];
            float y = (float)xyz[1];
            float z = (float)xyz[2];
            glTexCoord2f(s,t1);
            glVertex3f(x,y,z);
            xyz = _model.xyz(lat2,lon,0.0);
            x = (float)xyz[0];
            y = (float)xyz[1];
            z = (float)xyz[2];
            glTexCoord2f(s,t2);
            glVertex3f(x,y,z);
          }
          glEnd();
        }
        glEndList();
      }

      // Prepare to draw the texture.
      glShadeModel(GL_FLAT);
      glEnable(GL_TEXTURE_2D);
      glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);
      glBindTexture(GL_TEXTURE_2D,_textureName.name());
      glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

      // Draw the image.
      glPushMatrix();
      glScaled(_scale,_scale,_scale);
      glCallList(_displayList.list());
      glPopMatrix();

      // Done with the texture.
      glBindTexture(GL_TEXTURE_2D,0);
      glDisable(GL_TEXTURE_2D);
    }

    private OblateModel _model;
    private SampledImage _image;
    private double _scale;
    private GlDisplayList _displayList;
    private GlTextureName _textureName;
    private BoundingSphere _bs;
  }

  private static class GridNode extends Node {

    public GridNode(
      OblateModel model, Sampling slon, Sampling slat, double scale) 
    {
      _model = model;
      _slon = slon;
      _slat = slat;
      _scale = scale;
      _bs = new BoundingSphere(0.0,0.0,0.0,_model.getRadius());
    }

    protected BoundingSphere computeBoundingSphere(boolean finite) {
      return _bs;
    }

    protected void draw(DrawContext dc) {

      // If necessary, make display list and texture.
      if (_displayList==null) {
        _displayList = new GlDisplayList();
        glNewList(_displayList.list(),GL_COMPILE);
        int nlon = _slon.getCount();
        int nlat = _slat.getCount();
        float flon = (float)_slon.getFirst();
        float llon = (float)_slon.getLast();
        float flat = (float)_slat.getFirst();
        float llat = (float)_slat.getLast();
        float dlonl = 2.0f;
        int nlonl = 2+(int)((llon-flon)/dlonl);
        dlonl = (llon-flon)/nlonl;
        for (int ilat=0; ilat<nlat; ++ilat) {
          float lat = (float)_slat.getValue(ilat);
          glBegin(GL_LINE_STRIP);
          for (int i=0; i<=nlonl; ++i) {
            float lon = flon+i*dlonl;
            double[] xyz = _model.xyz(lat,lon,0.0);
            float x = (float)xyz[0];
            float y = (float)xyz[1];
            float z = (float)xyz[2];
            glVertex3f(x,y,z);
          }
          glEnd();
        }
        float dlatl = 2.0f;
        int nlatl = 2+(int)((llat-flat)/dlatl);
        dlatl = (llat-flat)/nlatl;
        for (int ilon=0; ilon<nlon; ++ilon) {
          float lon = (float)_slon.getValue(ilon);
          glBegin(GL_LINE_STRIP);
          for (int i=0; i<=nlatl; ++i) {
            float lat = flat+i*dlatl;
            double[] xyz = _model.xyz(lat,lon,0.0);
            float x = (float)xyz[0];
            float y = (float)xyz[1];
            float z = (float)xyz[2];
            glVertex3f(x,y,z);
          }
          glEnd();
        }
        glEndList();
      }

      // Draw the grid.
      glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
      glPushMatrix();
      glScaled(_scale,_scale,_scale);
      glCallList(_displayList.list());
      glPopMatrix();
    }

    private OblateModel _model;
    private Sampling _slon,_slat;
    private double _scale;
    private GlDisplayList _displayList;
    private BoundingSphere _bs;
  }
}
