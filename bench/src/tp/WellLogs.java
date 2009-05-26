package tp;

import java.io.*;
import java.util.*;
import edu.mines.jtk.util.Check;

/**
 * Well log processing for Teapot Dome.
 * 
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.05.21
 */
public class WellLogs {
  public void load(String dirName) {
    //System.out.println("WellLogs.load: "+dirName);
    File dir = new File(dirName);
    Check.argument(dir.isDirectory(),dirName+" is a directory");
    File[] files = dir.listFiles();
    for (File file:files) {
      if (file.isDirectory()) {
        load(file.getPath());
      } else if (file.isFile()) {
        WellLog log = WellLog.load(file);
        if (log!=null)
          _logs.add(log);
      }
    }
    System.out.println("WellLogs.load: nlogs = "+_logs.size());
    removeDuplicateLogs();
    System.out.println("WellLogs.load: nlogs = "+_logs.size());
  }
  private void removeDuplicateLogs() {
    HashSet<String> ids = new HashSet<String>();
    HashSet<Integer> rs = new HashSet<Integer>();
    for (int i=0; i<_logs.size(); ++i) {
      WellLog log = _logs.get(i);
      if (ids.contains(log.id))
        rs.add(i);
      else
        ids.add(log.id);
    }
    ArrayList<WellLog> logs = new ArrayList<WellLog>();
    for (int i=0; i<_logs.size(); ++i) {
      WellLog log = _logs.get(i);
      if (!rs.contains(i))
        logs.add(log);
    }
    _logs = logs;
  }
  private ArrayList<WellLog> _logs = new ArrayList<WellLog>();
}
