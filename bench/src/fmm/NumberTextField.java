/****************************************************************************
Copyright (c) 2008, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package fmm;

import java.awt.*;
import java.text.ParseException;
import javax.swing.*;
import javax.swing.text.*;

import edu.mines.jtk.util.StringUtil;

/**
 * A text field for numbers, either integer or floating-point values.
 * A number text field permits only numeric entries. Values may also be 
 * constrained by lower and upper bounds and to be integers.
 * <p>
 * The number displayed in a number text field is always a number or a
 * part of a number not yet completed. A number text field does not 
 * permit characters to be entered that could not be part of a number. 
 * <p>
 * When lower and upper bounds are specified, a number text field with 
 * keyboard focus will attempt to retain that focus until the displayed
 * value entered via the keyboard is within bounds.
 * <p>
 * Likewise, when constrained to be an integer, a number text field will
 * attempt to retain keyboard focus until the displayed value is an integer.
 * @author Dave Hale, Colorado School of Mines
 * @version 2008.06.21
 */
public class NumberTextField extends JFormattedTextField {
  private static final long serialVersionUID = 1L;

  /**
   * Constructs a number text field with no bounds.
   */
  public NumberTextField() {
    this(-Double.MAX_VALUE,Double.MAX_VALUE);
  }

  /**
   * Constructs a number text field with specified bounds on values.
   * @param vmin the minimum value.
   * @param vmax the maximum value.
   */
  public NumberTextField(double vmin, double vmax) {
    this(vmin,vmax,false);
  }

  /**
   * Constructs a number text field with specified constraints.
   * @param vmin the minimum value.
   * @param vmax the maximum value.
   * @param vint true, for values constrained to integers; false, otherwise.
   */
  public NumberTextField(double vmin, double vmax, final boolean vint) {
    super(new CustomFormatter(vmin,vmax,vint));
    _vmin = vmin;
    _vmax = vmax;
    _vint = vint;
    setInputVerifier(new InputVerifier() {
      public boolean verify(JComponent component) {
        String s = getText();
        double v;
        try {
          v = (_vint)?Integer.parseInt(s):Double.parseDouble(s);
        } catch (NumberFormatException e) {
          return false;
        }
        return _vmin<=v && v<=_vmax;
      }
    });
  }

  /**
   * Sets the printf-style format used to display the value in this field.
   * @param format the format.
   */
  public void setFormat(String format) {
    CustomFormatter cf = (CustomFormatter)getFormatter();
    cf.setFormat(format);
    repaint();
  }

  /**
   * Sets the min-max range of values. If necessary, this method modifies
   * the current value to be within the specified range of values.
   * @param vmin the minimum value.
   * @param vmax the maximum value.
   */
  public void setValueRange(double vmin, double vmax) {
    _vmin = vmin;
    _vmax = vmax;
    CustomFormatter cf = (CustomFormatter)getFormatter();
    cf.setValueRange(_vmin,_vmax);
    if (getDouble()<_vmin)
      setDouble(_vmin);
    if (getDouble()>_vmax)
      setDouble(_vmax);
  }

  /**
   * Sets the value of this number text field.
   * @param object the value object; must be a Number.
   */
  public void setValue(Object object) {
    Number n = (Number)object;
    double v = n.doubleValue();
    if (v<_vmin) v = _vmin;
    if (v>_vmax) v = _vmax;
    if (_vint)
      v = (int)v;
    super.setValue(v);
  }

  /**
   * Sets the value of this number text field as a double.
   * @param v the value.
   */
  public void setDouble(double v) {
    setValue(new Double(v));
  }

  /**
   * Sets the value of this number text field as a float.
   * @param v the value.
   */
  public void setFloat(float v) {
    setDouble(v);
  }

  /**
   * Sets the value of this number text field as an int.
   * @param v the value.
   */
  public void setInt(int v) {
    setDouble(v);
  }

  /**
   * Gets the value of this number text field as a double.
   * @return the value.
   */
  public double getDouble() {
    Number value = (Number)getValue();
    return value.doubleValue();
  }

  /**
   * Gets the value of this number text field as a float.
   * @return the value.
   */
  public float getFloat() {
    return (float)getDouble();
  }

  /**
   * Gets the value of this number text field as an int.
   * @return the value.
   */
  public int getInt() {
    return (int)getDouble();
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private double _vmin = -Double.MAX_VALUE; // min value
  private double _vmax = Double.MAX_VALUE; // max value
  private boolean _vint = false; // true for only integer values

  private static double constrainValue(
    double v, double vmin, double vmax, boolean vint) 
  {
    if (v<vmin) v = vmin;
    if (v>vmax) v = vmax;
    if (vint)
      v = (int)v;
    return v;
  }

  // Formatter that uses a specified printf-style format for display.
  // Also removes any insignificant trailing zeros or decimal point.
  private static class CustomFormatter extends DefaultFormatter {
    public CustomFormatter(double vmin, double vmax, boolean vint) {
      this(vmin,vmax,vint,"%1.6g");
    }
    public CustomFormatter(
      double vmin, double vmax, boolean vint, String format) 
    {
      _vmin = vmin;
      _vmax = vmax;
      _vint = vint;
      _format = format;
      _filter = new CustomFilter(vint);
      setOverwriteMode(true);
      setAllowsInvalid(false);
    }
    public void setFormat(String format) {
      _format = format;
    }
    public void setValueRange(double vmin, double vmax) {
      _vmin = vmin;
      _vmax = vmax;
    }
    public String valueToString(Object v) throws ParseException {
      if (!(v instanceof Double)) 
        throw new ParseException("value is not a double",0);
      String s = String.format(_format,(Double)v);
      s = StringUtil.removeTrailingZeros(s);
      return s;
    }
    public Object stringToValue(String s) throws ParseException {
      Double v = null;
      try {
        v = (_vint)?Integer.parseInt(s):Double.parseDouble(s);
        if (v<_vmin || v>_vmax) throw new ParseException("out of range",0);
      } catch (NumberFormatException e) {
        throw new ParseException("cannot convert string to double",0);
      }
      return v;
    }
    protected DocumentFilter getDocumentFilter() {
      return _filter;
    }
    private CustomFilter _filter;
    private String _format;
    private double _vmin,_vmax;
    private boolean _vint;
  }

  // Ensures that the text in the field is always a valid part of a 
  // number. The trick here is to append a zero to the changed text.
  // (The changed text is the current text with the proposed change.)
  // If the changed text with zero appended is a valid number, then 
  // we permit the proposed change. Otherwise, we beep.
  private static class CustomFilter extends DocumentFilter {
    public CustomFilter(boolean vint) {
      _vint = vint;
    }
    public void insertString(
      FilterBypass fb, int off, String string, AttributeSet as)
      throws BadLocationException
    {
      String sc = getCurrentText(fb);
      StringBuilder sb = new StringBuilder(sc);
      String sn = sb.insert(off,string).toString();
      if (isPartOfNumber(sn)) {
        super.insertString(fb,off,string,as);
      } else {
        Toolkit.getDefaultToolkit().beep();
      }
    }
    public void replace(
      FilterBypass fb, int off, int len, String string, AttributeSet as)
      throws BadLocationException
    {
      String sc = getCurrentText(fb);
      StringBuilder sb = new StringBuilder(sc);
      String sn = sb.replace(off,off+len,string).toString();
      if (isPartOfNumber(sn)) {
        super.replace(fb,off,len,string,as);
      } else {
        Toolkit.getDefaultToolkit().beep();
      }
    }
    private boolean _vint;
    private String getCurrentText(FilterBypass fb) {
      Document d = fb.getDocument();
      String s = null;
      try {
        s = d.getText(0,d.getLength());
      } catch (BadLocationException e) {
        assert false:"exception not possible: "+e;
      }
      return s;
    }
    private boolean isPartOfNumber(String s) {
      s = s+"0";
      try {
        if (_vint) {
          int i = Integer.parseInt(s);
        } else {
          double d = Double.parseDouble(s);
        }
        return true;
      } catch (NumberFormatException e) {
        return false;
      }
    }
  }
}
