package scaniter;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Collections;
import java.util.StringTokenizer;

import modi.Constants;
import modi.Peak;
import modi.Spectrum;

public class MSMScan {
  static int minPeaksCount = 4;
  static double minMW = 8 * 57 + 18;

  private String title;
  private int scanNo;
  private double pmz;
  private double neutralMW;
  private int charge;
  private long offset;
  private Spectrum peaklist;
  private static double tolerance = Constants.massToleranceForDenovo;

  public MSMScan(String title, double pmz, int charge) {
    this.scanNo = Integer.valueOf(title);
    this.pmz = pmz;
    this.charge = charge;
    this.neutralMW = (pmz - Constants.Proton) * charge;
  }

  public MSMScan(String title, int sn, double pmz, int charge) {
    this.title = title;
    this.scanNo = sn;
    this.pmz = pmz;
    this.charge = charge;
    this.neutralMW = (pmz - Constants.Proton) * charge;
  }


  public String getTitle() {
    return title;
  }

  public int getScanNumber() {
    return scanNo;
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

  public int getPeaksCount() {
    return peaklist.size();
  }

  public long getOffset() {
    return offset;
  }

  public Spectrum getSpectrum() {
    return peaklist;
  }

  public void setOffset(long offset) {
    this.offset = offset;
  }

  public void readPeakList(RandomAccessFile in) throws IOException {

    if (Constants.rangeForIsotopeIncrement != 0)
      Constants.NoOfC13 =
          (int) Math.ceil(neutralMW / Constants.rangeForIsotopeIncrement);
    if (Constants.PPMTolerance != 0)
      Constants.precursorAccuracy =
          Constants.PPMtoDalton(neutralMW, Constants.PPMTolerance);
    Constants.precursorTolerance =
        Constants.precursorAccuracy + Constants.NoOfC13 * Constants.IsotopeSpace;

    String s;
    ArrayList<RawP> rawPL = new ArrayList<RawP>();

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

    if (Constants.iTRAQSearch) {

      double reporterIon114 = 114.1105;
      double reporterIon115 = 115.1074;
      double reporterIon116 = 116.1107;
      double reporterIon117 = 117.1141;
      double reporterIonTolerance = 0.1;

      double isobaricTag = 144.102063;
      double isobarictagTolerance = 0.005;

      double leftTagMass = isobaricTag + Constants.Proton;
      double
          oneChargeRightTagMass =
          this.neutralMW - isobaricTag + Constants.Proton;
      double
          twoChargeRightTagMass =
          (this.neutralMW - isobaricTag + (Constants.Proton * 2)) / 2;

      double isotopeDelta = 1.0034 / this.charge;
      int isotopeRatio = 5;

      for (int j = 0; j < rawPL.size(); j++) {
        double ionMZ = rawPL.get(j).mz;
        double ionIT = rawPL.get(j).it;

        if (Math.abs(ionMZ - reporterIon114) < reporterIonTolerance || Math.abs(ionMZ - reporterIon115) < reporterIonTolerance || Math.abs(ionMZ - reporterIon116) < reporterIonTolerance || Math.abs(ionMZ - reporterIon117) < reporterIonTolerance) {
          rawPL.remove(j);
          j--;
        }
        else if ((Math.abs(ionMZ - leftTagMass) < isobarictagTolerance) || (this.charge == 2 && Math.abs(ionMZ - oneChargeRightTagMass) < isobarictagTolerance)) {
          rawPL.remove(j);
          int ISO = j;
          double prevISO = ionIT;
          boolean isoDecent = false;
          double targetmz = ionMZ + isotopeDelta;
          for (int nth_iso = 1; ; nth_iso++) {
            int H = -1;
            double imax = 0;
            for (int i = ISO; i < rawPL.size(); i++) {
              if (rawPL.get(i).mz > targetmz + isobarictagTolerance) {
                break;
              }
              else if (Math.abs(rawPL.get(i).mz - targetmz) < isobarictagTolerance) {
                if (rawPL.get(i).it > imax) {
                  imax = rawPL.get(i).it;
                  H = i;
                }
              }
            }

            if (H == -1) {
              break;
            }
            if (prevISO < imax / isotopeRatio) {
              break;
            }
            if (isoDecent && (prevISO < imax)) {
              break;
            }
            if (prevISO > imax) {
              isoDecent = true;
            }

            targetmz = rawPL.get(H).mz;
            prevISO = imax;
            rawPL.remove(H);
            ISO = H;
          }
          j--;
        }
        else if (this.charge == 3 && Math.abs(ionMZ - twoChargeRightTagMass) < isobarictagTolerance) {
          double okmax = 0;
          double targetmz = ionMZ - isotopeDelta;
          for (int i = j - 1; i > -1; i--) {
            if (rawPL.get(i).mz < targetmz - isobarictagTolerance)
              break;
            else if (Math.abs(rawPL.get(i).mz - targetmz) < isobarictagTolerance) {
              if (rawPL.get(i).it > okmax) {
                okmax = rawPL.get(i).it;
              }
            }
          }
          if (okmax > ionIT / isotopeRatio) {
            continue;
          }

          rawPL.remove(j);
          int ISO = j;
          double prevISO = ionIT;
          boolean isoDecent = false;
          targetmz = ionMZ + isotopeDelta;
          for (int nth_iso = 1; ; nth_iso++) {
            int H = -1;
            double imax = 0;
            for (int i = ISO; i < rawPL.size(); i++) {
              if (rawPL.get(i).mz > targetmz + isobarictagTolerance) {
                break;
              }
              else if (Math.abs(rawPL.get(i).mz - targetmz) < isobarictagTolerance) {
                if (rawPL.get(i).it > imax) {
                  imax = rawPL.get(i).it;
                  H = i;
                }
              }
            }

            if (H == -1) {
              break;
            }
            if (prevISO < imax / isotopeRatio) {
              break;
            }
            if (isoDecent && (prevISO < imax)) {
              break;
            }
            if (prevISO > imax) {
              isoDecent = true;
            }

            targetmz = rawPL.get(H).mz;
            prevISO = imax;
            rawPL.remove(H);
            ISO = H;
          }
          j--;
        }
      }

    }

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

    if (spectrum.size() < minPeaksCount || neutralMW < minMW)
      peaklist = null;
    else
      peaklist = spectrum;
  }

  public void readPeakList(double[][] pList) throws IOException {

    if (Constants.rangeForIsotopeIncrement != 0)
      Constants.NoOfC13 =
          (int) Math.ceil(neutralMW / Constants.rangeForIsotopeIncrement);
    if (Constants.PPMTolerance != 0)
      Constants.precursorAccuracy =
          Constants.PPMtoDalton(neutralMW, Constants.PPMTolerance);
    Constants.precursorTolerance =
        Constants.precursorAccuracy + Constants.NoOfC13 * Constants.IsotopeSpace;

    ArrayList<RawP> rawPL = new ArrayList<RawP>();
    for (int j = 0; j < pList[0].length; j++) {
      rawPL.add(new RawP(pList[0][j], pList[1][j]));
    }
    Collections.sort(rawPL);

    if (Constants.iTRAQSearch) {

      double reporterIon114 = 114.1105;
      double reporterIon115 = 115.1074;
      double reporterIon116 = 116.1107;
      double reporterIon117 = 117.1141;
      double reporterIonTolerance = 0.1;

      double isobaricTag = 144.102063;
      double isobarictagTolerance = 0.005;

      double leftTagMass = isobaricTag + Constants.Proton;
      double
          oneChargeRightTagMass =
          this.neutralMW - isobaricTag + Constants.Proton;
      double
          twoChargeRightTagMass =
          (this.neutralMW - isobaricTag + (Constants.Proton * 2)) / 2;

      double isotopeDelta = 1.0034 / this.charge;
      int isotopeRatio = 5;

      for (int j = 0; j < rawPL.size(); j++) {
        double ionMZ = rawPL.get(j).mz;
        double ionIT = rawPL.get(j).it;

        if (Math.abs(ionMZ - reporterIon114) < reporterIonTolerance || Math.abs(ionMZ - reporterIon115) < reporterIonTolerance || Math.abs(ionMZ - reporterIon116) < reporterIonTolerance || Math.abs(ionMZ - reporterIon117) < reporterIonTolerance) {
          rawPL.remove(j);
          j--;
        }
        else if ((Math.abs(ionMZ - leftTagMass) < isobarictagTolerance) || (this.charge == 2 && Math.abs(ionMZ - oneChargeRightTagMass) < isobarictagTolerance)) {
          rawPL.remove(j);
          int ISO = j;
          double prevISO = ionIT;
          boolean isoDecent = false;
          double targetmz = ionMZ + isotopeDelta;
          for (int nth_iso = 1; ; nth_iso++) {
            int H = -1;
            double imax = 0;
            for (int i = ISO; i < rawPL.size(); i++) {
              if (rawPL.get(i).mz > targetmz + isobarictagTolerance) {
                break;
              }
              else if (Math.abs(rawPL.get(i).mz - targetmz) < isobarictagTolerance) {
                if (rawPL.get(i).it > imax) {
                  imax = rawPL.get(i).it;
                  H = i;
                }
              }
            }

            if (H == -1) {
              break;
            }
            if (prevISO < imax / isotopeRatio) {
              break;
            }
            if (isoDecent && (prevISO < imax)) {
              break;
            }
            if (prevISO > imax) {
              isoDecent = true;
            }

            targetmz = rawPL.get(H).mz;
            prevISO = imax;
            rawPL.remove(H);
            ISO = H;
          }
          j--;
        }
        else if (this.charge == 3 && Math.abs(ionMZ - twoChargeRightTagMass) < isobarictagTolerance) {
          double okmax = 0;
          double targetmz = ionMZ - isotopeDelta;
          for (int i = j - 1; i > -1; i--) {
            if (rawPL.get(i).mz < targetmz - isobarictagTolerance)
              break;
            else if (Math.abs(rawPL.get(i).mz - targetmz) < isobarictagTolerance) {
              if (rawPL.get(i).it > okmax) {
                okmax = rawPL.get(i).it;
              }
            }
          }
          if (okmax > ionIT / isotopeRatio) {
            continue;
          }

          rawPL.remove(j);
          int ISO = j;
          double prevISO = ionIT;
          boolean isoDecent = false;
          targetmz = ionMZ + isotopeDelta;
          for (int nth_iso = 1; ; nth_iso++) {
            int H = -1;
            double imax = 0;
            for (int i = ISO; i < rawPL.size(); i++) {
              if (rawPL.get(i).mz > targetmz + isobarictagTolerance) {
                break;
              }
              else if (Math.abs(rawPL.get(i).mz - targetmz) < isobarictagTolerance) {
                if (rawPL.get(i).it > imax) {
                  imax = rawPL.get(i).it;
                  H = i;
                }
              }
            }

            if (H == -1) {
              break;
            }
            if (prevISO < imax / isotopeRatio) {
              break;
            }
            if (isoDecent && (prevISO < imax)) {
              break;
            }
            if (prevISO > imax) {
              isoDecent = true;
            }

            targetmz = rawPL.get(H).mz;
            prevISO = imax;
            rawPL.remove(H);
            ISO = H;
          }
          j--;
        }
      }

    }

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

    if (spectrum.size() < minPeaksCount || neutralMW < minMW)
      peaklist = null;
    else
      peaklist = spectrum;
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
}
