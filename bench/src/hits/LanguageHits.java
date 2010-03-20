package hits;

import java.io.*;
import java.net.*;
import java.util.Date;

/**
 * Statistics for programming languages based on search engine hit counts.
 * For each programming language lang, this program uses the search phrase 
 * "lang programming" "field".
 * The first part of this phrase is that used by TIOBE (www.tiobe.com).
 * The second part restricts the search to a specific field, such as
 * "geophysics" or "high performance computing".
 * <p>
 * Note that the GSERP method violates Google's Terms of Service, as
 * this method accesses Google's search service through an automated
 * means that is not specifically allowed. The GOOGLE method uses the
 * web API provided by Google, which is specifically allowed. 
 * <p> 
 * Likewise, the YAHOO and BING methods use Yahoo's and Microsoft's web 
 * APIs for their respective search engines. These two methods also 
 * require developer IDs (API keys) that are freely available, and these
 * IDs must be specified as command-line arguments: args[0] for Yahoo,
 * args[1] for Bing.
 * @author Dave Hale, Colorado School of Mines
 * @version 2010.02.26
 */
public class LanguageHits {
  private enum Engine { GSERP, GOOGLE, YAHOO, BING };
  private static final String[] LANGS = {
    "Java","C","C++","Python","MATLAB","Fortran",
  };
  private static final String[] FIELDS = {
    "","scientific","engineering",
    "geophysics", "physics","geology",
    "petroleum engineering","electrical engineering",
    "high performance computing",
    "chemistry",
  };
  private static final int NLANG = LANGS.length;
  private static final int NFIELD = FIELDS.length;
  private static String YAHOO_API_KEY; // get from args[0]
  private static String BING_API_KEY;  // get from args[1]

  public static void main (String[] args) throws IOException{
    if (args.length>0) YAHOO_API_KEY = args[0];
    if (args.length>1)  BING_API_KEY = args[1];
    System.out.println(new Date());
    Engine[] engines = {
      Engine.GSERP, 
      Engine.GOOGLE,
      Engine.YAHOO,
      Engine.BING};
    for (Engine engine:engines) {
      if (engine==Engine.YAHOO && YAHOO_API_KEY==null) continue;
      if (engine==Engine.BING  &&  BING_API_KEY==null) continue;
      System.out.println(engine);
      int[][] c = countHits(engine);
      printTable(c);
    }
  }

  private static void printTable(int[][] c) {
    int[] t = totalHits(c);
    double[][] f = fracHits(c);
    System.out.println("Counts");
    printHeadings();
    for (int ilang=0; ilang<NLANG; ++ilang) {
      System.out.printf("%7s",LANGS[ilang]);
      for (int ifield=0; ifield<NFIELD; ++ifield) {
        System.out.print(","+c[ilang][ifield]);
      }
      System.out.println();
    }
    System.out.print("Total");
    for (int ifield=0; ifield<NFIELD; ++ifield)
      System.out.print(","+t[ifield]);
    System.out.println();
    System.out.println("Percentages");
    printHeadings();
    for (int ilang=0; ilang<NLANG; ++ilang) {
      System.out.printf("%7s",LANGS[ilang]);
      for (int ifield=0; ifield<NFIELD; ++ifield) {
        double p = 100.0*f[ilang][ifield];
        System.out.printf(",%5.1f",p);
      }
      System.out.println();
    }
  }
  private static void printHeadings() {
    System.out.print("language");
    for (String field:FIELDS)
      System.out.print(","+field);
    System.out.println();
  }

  private static void sleep(int seconds) {
    try {
      Thread.currentThread().sleep(seconds*1000);
    } catch (Exception e) {
    }
  }
  private static int[][] countHits(Engine engine) {
    int[][] c = new int[NLANG][NFIELD];
    for (int ilang=0; ilang<NLANG; ++ilang) {
      for (int ifield=0; ifield<NFIELD; ++ifield) {
        if (engine==Engine.GOOGLE) {
          c[ilang][ifield] = countHitsGoogle(LANGS[ilang],FIELDS[ifield]);
        } else if (engine==Engine.YAHOO) {
          c[ilang][ifield] = countHitsYahoo(LANGS[ilang],FIELDS[ifield]);
        } else if (engine==Engine.BING) {
          c[ilang][ifield] = countHitsBing(LANGS[ilang],FIELDS[ifield]);
        } else if (engine==Engine.GSERP) {
          c[ilang][ifield] = countHitsGserp(LANGS[ilang],FIELDS[ifield]);
        }
        //sleep(2);
      }
    }
    return c;
  }
  private static int[] totalHits(int[][] c) {
    int[] t = new int[NFIELD];
    for (int ilang=0; ilang<NLANG; ++ilang)
      for (int ifield=0; ifield<NFIELD; ++ifield)
        t[ifield] += c[ilang][ifield];
    return t;
  }
  private static double[][] fracHits(int[][] c) {
    int[] t = totalHits(c);
    double[][] f = new double[NLANG][NFIELD];
    for (int ilang=0; ilang<NLANG; ++ilang)
      for (int ifield=0; ifield<NFIELD; ++ifield)
        f[ilang][ifield] = (double)c[ilang][ifield]/(double)t[ifield];
    return f;
  }
  private static String getPageText(URL url) {
    StringBuilder sb = new StringBuilder();
    try {
      HttpURLConnection uc = (HttpURLConnection)url.openConnection();
      uc.setRequestProperty("User-agent","Mozilla/5.0");
      /*
      int code = uc.getResponseCode();
      String message = uc.getResponseMessage();
      System.out.println("code="+code);
      System.out.println("message="+message);
      */
      InputStream is = new BufferedInputStream(uc.getInputStream());
      Reader r = new InputStreamReader(is);
      int c;
      while ((c=r.read())!=-1)
        sb.append((char)c);
    } catch(IOException e) {
      throw new RuntimeException(e);
    }
    return sb.toString();
  }

  ///////////////////////////////////////////////////////////////////////////
  // Google SERP (scrapes the Search Engine Results Page)
  private static int countHitsGserp(String lang, String field) {
    String COUNT_PREFIX = "> of about <b>";
    String COUNT_SUFFIX = "</b> for";
    URL url = makeUrlGserp(lang,field);
    String pt = getPageText(url);
    int i = pt.indexOf(COUNT_PREFIX)+COUNT_PREFIX.length();
    int j = pt.indexOf(COUNT_SUFFIX,i);
    if (i<0 || j<=i || j-i>20) 
      throw new RuntimeException("cannot grok the search results");
    String cs = pt.substring(i,j);
    return Integer.parseInt(cs.replaceAll(",",""));
  }
  private static URL makeUrlGserp(String lang, String field) {
    StringBuilder sb = new StringBuilder("http://www.google.com");
    sb.append("/search?hl=en&q=%22");
    sb.append(lang.replace("+","%2B").replace("#","%23"));
    sb.append("+programming%22+%22");
    sb.append(field.replace(" ","+"));
    sb.append("%22&aq=f&aqi=&aql=&oq=");
    URL url = null;
    try {
      url = new URL(sb.toString());
    } catch (MalformedURLException e) {
      throw new RuntimeException(e);
    }
    return url;
  }

  ///////////////////////////////////////////////////////////////////////////
  // Google
  private static int countHitsGoogle(String lang, String field) {
    String COUNT_PREFIX = "\"estimatedResultCount\":\"";
    String COUNT_SUFFIX = "\",";
    URL url = makeUrlGoogle(lang,field);
    String pt = getPageText(url);
    int i = pt.indexOf(COUNT_PREFIX)+COUNT_PREFIX.length();
    int j = pt.indexOf(COUNT_SUFFIX,i);
    if (i<0 || j<=i || j-i>20) 
      throw new RuntimeException("cannot grok the search results");
    String cs = pt.substring(i,j);
    return Integer.parseInt(cs.replaceAll(",",""));
  }
  private static URL makeUrlGoogle(String lang, String field) {
    StringBuilder sb = new StringBuilder("http://ajax.googleapis.com");
    sb.append("/ajax/services/search/web?v=1.0&q=%22");
    sb.append(lang.replace("+","%2B").replace("#","%23"));
    sb.append("+programming%22+%22");
    sb.append(field.replace(" ","+"));
    sb.append("%22");
    URL url = null;
    try {
      url = new URL(sb.toString());
    } catch (MalformedURLException e) {
      throw new RuntimeException(e);
    }
    return url;
  }

  ///////////////////////////////////////////////////////////////////////////
  // Yahoo
  private static int countHitsYahoo(String lang, String field) {
    String COUNT_PREFIX = "\"totalhits\":\"";
    String COUNT_SUFFIX = "\",";
    URL url = makeUrlYahoo(lang,field);
    //System.out.println(url);
    String pt = getPageText(url);
    int i = pt.indexOf(COUNT_PREFIX)+COUNT_PREFIX.length();
    int j = pt.indexOf(COUNT_SUFFIX,i);
    if (i<0 || j<=i || j-i>20) 
      throw new RuntimeException("cannot grok the search results");
    String cs = pt.substring(i,j);
    return Integer.parseInt(cs.replaceAll(",",""));
  }
  private static URL makeUrlYahoo(String lang, String field) {
    StringBuilder sb = new StringBuilder("http://boss.yahooapis.com");
    sb.append("/ysearch/web/v1/%22");
    sb.append(lang.replace("+","%2B").replace("#","%23"));
    sb.append("+programming%22+%22");
    sb.append(field.replace(" ","+"));
    sb.append("%22?appid=");
    sb.append(YAHOO_API_KEY);
    sb.append("&format=json");
    URL url = null;
    try {
      url = new URL(sb.toString());
    } catch (MalformedURLException e) {
      throw new RuntimeException(e);
    }
    return url;
  }

  ///////////////////////////////////////////////////////////////////////////
  // Bing
  private static int countHitsBing(String lang, String field) {
    String COUNT_PREFIX = "{\"Total\":";
    String COUNT_SUFFIX = ",";
    URL url = makeUrlBing(lang,field);
    //System.out.println(url);
    String pt = getPageText(url);
    int i = pt.indexOf(COUNT_PREFIX)+COUNT_PREFIX.length();
    int j = pt.indexOf(COUNT_SUFFIX,i);
    if (i<0 || j<=i || j-i>20) 
      throw new RuntimeException("cannot grok the search results");
    String cs = pt.substring(i,j);
    return Integer.parseInt(cs.replaceAll(",",""));
  }
  private static URL makeUrlBing(String lang, String field) {
    StringBuilder sb = new StringBuilder("http://api.search.live.net");
    sb.append("/json.aspx?AppId=");
    sb.append(BING_API_KEY);
    sb.append("&Query=%22");
    sb.append(lang.replace("+","%2B").replace("#","%23"));
    sb.append("+programming%22+%22");
    sb.append(field.replace(" ","+"));
    sb.append("&Sources=Web");
    URL url = null;
    try {
      url = new URL(sb.toString());
    } catch (MalformedURLException e) {
      throw new RuntimeException(e);
    }
    return url;
  }
}
