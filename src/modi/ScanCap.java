package modi;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Collections;
import java.util.StringTokenizer;

public class ScanCap implements Comparable<ScanCap> {
  private String title;
  private double pmz;
  private double neutralMW;
  private int charge;
  private long offset;

  private static double tolerance = Constants.massToleranceForDenovo;

  public ScanCap(String title, double pmz, int charge) {
    this.title = title;
    this.pmz = pmz;
    this.charge = charge;
    this.neutralMW = (pmz - Constants.Proton) * charge;
  }

  public void setOffset(long offset) {
    this.offset = offset;
  }

  public String getTitle() {
    return title;
  }

  public double getObservedMW() {
    return neutralMW;
  }

  public double getPMZ() {
    return pmz;
  }

  public int getCharge() {
    return charge;
  }

  public long getOffset() {
    return offset;
  }

  public Spectrum getSpectrum(RandomAccessFile in) throws IOException {

    if (Constants.rangeForIsotopeIncrement != 0) {
      Constants.NoOfC13 =
          (int) Math.ceil(neutralMW / Constants.rangeForIsotopeIncrement);
    }

    if (Constants.PPMTolerance != 0) {
      Constants.precursorAccuracy =
          Constants.PPMtoDalton(neutralMW, Constants.PPMTolerance);
    }
    Constants.precursorTolerance =
        Constants.precursorAccuracy + Constants.NoOfC13 * Constants.IsotopeSpace;

    String s;

    ArrayList<RawP> rawPL = new ArrayList<RawP>();
    in.seek(this.offset);
    while ((s = in.readLine()) != null) {
      StringTokenizer token = new StringTokenizer(s);
      if (token.countTokens() > 1) {
        if (!Character.isDigit(s.charAt(0)))
          break;
        rawPL.add(new RawP(Double.parseDouble(token.nextToken()), Double.parseDouble(token.nextToken())));
      }
      else
        break;
    }
    Collections.sort(rawPL);

    int index = 0;
    Spectrum spectrum = new Spectrum(this.pmz, this.charge, this.title);

    double basePeakIntensity = 0, TIC = 0;
    double tarMass = 0, tarInten = 0;
    for (RawP rp : rawPL) {
      double mass = rp.mz;
      double intensity = rp.it;

      if (intensity <= 0 || mass <= 0)
        continue;
      if (mass > neutralMW)
        continue;
      if (Math.abs(mass - pmz) < 5.)
        continue;

      if ((mass - tarMass) < tolerance) {
        double sum = tarInten + intensity;
        tarMass = tarMass * (tarInten / sum) + mass * (intensity / sum);
        tarInten += intensity;
        spectrum.get(index - 1).set(tarMass, tarInten);
      }
      else {
        spectrum.add(new Peak(index++, mass, intensity));
        tarMass = mass;
        tarInten = intensity;
      }
      TIC += intensity;
      if (tarInten > basePeakIntensity)
        basePeakIntensity = tarInten;
    }
    spectrum.setExtraInformation(basePeakIntensity, TIC);

    Constants.gapTolerance = Constants.fragmentTolerance * 2;
    Constants.nonModifiedDelta =
        (Constants.precursorTolerance < Constants.massToleranceForDenovo) ? Constants.precursorTolerance : Constants.massToleranceForDenovo;

    if (Constants.precursorTolerance > Constants.gapTolerance)
      Constants.gapTolerance = Constants.precursorTolerance;

    return spectrum;
  }

  private class RawP implements Comparable<RawP> {
    double mz;
    double it;

    public RawP(double m, double i) {
      mz = m;
      it = i;
    }

    public int compareTo(RawP p) {
      if (mz > p.mz)
        return 1;
      else if (mz < p.mz)
        return -1;
      else
        return 0;
    }
  }

  public int compareTo(ScanCap s) {
    if (this.neutralMW > s.neutralMW)
      return 1;
    else if (this.neutralMW < s.neutralMW)
      return -1;

    if (this.charge > s.charge)
      return 1;
    else if (this.charge < s.charge)
      return -1;
    else
      return 0;
  }
}
