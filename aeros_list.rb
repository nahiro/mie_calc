#!/usr/bin/ruby

LIST = 'list'
WREF = 550.00
WNEX = 514.50
EPSILON = 1.0e-14

File.readlines(LIST).each{|line|
  flag = 0
  begin
    item = line.chomp.split
    dnam = item[0]
    ncol = item.size
    if(!File.directory?(dnam))
      break
    end
    if(/^(\d\d)(\d\d)(\d\d)/ !~ dnam)
      $stderr.puts "Invalid DATE >>> #{dnam}"
      break
    end
    yr = $1.to_i
    mo = $2.to_i
    dy = $3.to_i
    if(yr<7 || yr>9 || mo<1 || mo>12 || dy<1 || dy>31)
      $stderr.puts "Invalid DATE >>> yr=#{yr},mo=#{mo},dy=#{dy}"
      break
    end
    torg = Time.local(yr,mo,dy)
    wval = 0.85170
    if(mo>=5 && mo<=10)
      wval = 2.92231
    end
    ino = -1
    rno = ""
    dir = ""
    pnam = ""
    onam = ""
    inam = ""
    99.downto(1){|i|
      rno_tmp = sprintf("%02d",i)
      dir_tmp = "#{dnam}/#{rno_tmp}"
      if(File.directory?(dir_tmp))
        pnam_tmp = "#{dir_tmp}/mie_pval.dat"
        onam_tmp = "#{dir_tmp}/mie_out1.dat"
        inam_tmp = "#{dir_tmp}/input.dat"
        if(File.exist?(pnam_tmp) && File.exist?(onam_tmp) && File.exist?(inam_tmp))
          ino = i
          rno = rno_tmp
          dir = dir_tmp
          onam = onam_tmp
          pnam = pnam_tmp
          inam = inam_tmp
          break
        end
      end
    }
    if(ino < 0)
      break
    end
    # Read mie_pval.dat
    lines = File.readlines(pnam)
    if(lines.size <= 10)
      $stderr.puts "Error, lines.size=#{lines.size} (#{pnam})"
      break
    end
    err = 0
    pval = []
    0.upto(10){|j|
      temp = lines[j]
      item = temp.chomp.split
      if(item.size != 2)
        $stderr.puts "Error, item.size=#{item.size} (#{pnam})"
        err = 1
        break
      end
      if(item[0].to_i != j)
        $stderr.puts "Error, item[0]=#{item[0]}, j=#{j} (#{pnam})"
        err = 1
        break
      end
      pval << item[1].to_f
    }
    if(err != 0)
      break
    end
    vmie = pval[0]
    wscl = pval[1]
    ncmp = pval[2]
    lef1 = pval[3]
    lsg1 = pval[4]
    mxr2 = pval[5]
    lef2 = pval[6]
    lsg2 = pval[7]
    mxr3 = pval[8]
    lef3 = pval[9]
    lsg3 = pval[10]
    wcm1 = 1.0/(1.0+1.0/10.0**mxr2+1.0/10.0**mxr3)
    wcm2 = wcm1/10.0**mxr2
    wcm3 = wcm1/10.0**mxr3
    # Read mie_out1.dat
    lines = File.readlines(onam)
    if(lines.size < 1)
      $stderr.puts "Error, lines.size=#{lines.size} (#{onam})"
      break
    end
    err = 0
    cext = 0.0
    omeg = 0.0
    asym = 0.0
    lines.each{|temp|
      item = temp.chomp.split
      if(item.size != 5)
        $stderr.puts "Error, item.size=#{item.size} (#{onam})"
        err = 1
        break
      end
      w = item[0].to_f
      e = item[2].to_f
      o = item[3].to_f
      a = item[4].to_f
      if((w-WNEX).abs < EPSILON)
        cext = e
      end
      if((w-WREF).abs < EPSILON)
        if((e-1.0).abs > EPSILON)
          $stderr.puts "Error, w=#{w}, e=#{e} (#{onam})"
          err = 1
          break
        end
        omeg = o
        asym = a
      end
    }
    if(err != 0)
      break
    end
    pang = -Math.log(cext)/Math.log(WNEX/WREF)
    # Read input.dat
    lines = File.readlines(inam)
    if(lines.size < 1)
      $stderr.puts "Error, lines.size=#{lines.size} (#{inam})"
      break
    end
    err = 0
    tdif = 0.0
    lines.each{|temp|
      item = temp.chomp.split
      if(item.size != 19)
        $stderr.puts "Error, item.size=#{item.size} (#{inam})"
        err = 1
        break
      end
      yr = item[0].to_i
      mo = item[1].to_i
      dy = item[2].to_i
      hr = item[3].to_i
      mi = item[4].to_i
      sc = item[5].to_i
      tdif += (Time.local(yr,mo,dy,hr,mi,sc)-torg)
    }
    if(err != 0)
      break
    end
    tdif /= lines.size
    tavg = torg+tdif
    # Output results
    fp = $stdout
    if(ncol != 1)
      fp = $stderr
    else
      flag = 1
    end
    fp.printf("%-10s %4.1f %s %5.2f %5.3f %5.3f   %6.1f %6.3f   %7.5f %7.5f %7.5f   %7.4f %7.4f %7.4f\n",
               dnam,tavg.hour+tavg.min/60.0+tavg.sec/3600.0,rno,pang,omeg,asym,
               vmie,wval*wscl,wcm1,wcm2,wcm3,10.0**lef1,10.0**lef2,10.0**lef3)
  end while(nil)
  if(flag == 0)
    puts line
  end
}
