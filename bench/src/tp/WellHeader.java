package tp;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;

/**
 * A well header from Teapot Dome.
 * 
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.05.26
 */
public class WellHeader {
  public long id; // unique 12-digit API well number
  public double xe,yn,ze; // easting, northing, and elevation (ft)

  /**
   * Collection of well header data.
   */
  public static class Data {

    /**
     * Constructs well header data from the specified file.
     * The file is in Excel CSV (comma-delimited) format.
     * @param fileName file containing well header data.
     */
    public Data(String fileName) {
      try {
        FileInputStream fis = new FileInputStream(fileName);
        Scanner s = new Scanner(fis);
        while (s.hasNextLine()) {
          String line = s.nextLine();
          String[] fields = line.split(",");
          if (fields.length<10)
            continue;
          long id = WellLog.idFromString(fields[0]);
          if (id<0)
            continue;
          WellHeader wh = new WellHeader();
          wh.id = id;
          try {
            wh.xe = Double.parseDouble(fields[5]);
            wh.yn = Double.parseDouble(fields[4]);
            wh.ze = Double.parseDouble(fields[9]);
            add(wh);
          } catch (NumberFormatException e) {
            // do nothing if well header is missing something
          }
        }
        s.close();
      } catch (IOException e) {
        throw new RuntimeException(e);
      }
    }

    /**
     * Adds the specified well header.
     * @param wh the well header.
     */
    public void add(WellHeader wh) {
      _data.put(wh.id,wh);
    }

    /**
     * Returns the number of well headers.
     * @return the number of well headers.
     */
    public int size() {
      return _data.size();
    }

    /**
     * Gets the well header for the specified well id.
     * @param id the well id.
     * @return the well header; null, if none.
     */
    public WellHeader get(long id) {
      return _data.get(id);
    }

    /**
     * Gets all well headers.
     * @return list of well headers.
     */
    public List<WellHeader> getAll() {
      return new ArrayList<WellHeader>(_data.values());
    }

    /**
     * Prints summary information for these well header data.
     */
    public void printInfo() {
      for (WellHeader wh:getAll())
        System.out.println("id="+wh.id+" xe="+wh.xe+" yn="+wh.yn+" ze="+wh.ze);
      System.out.println("number of headers = "+size());
    }
    private Map<Long,WellHeader> _data = new HashMap<Long,WellHeader>();
  }
}
