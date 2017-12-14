#!/usr/bin/ruby
CALDIR = "#{ENV['HOME']}/work/MS720/data/cal"
DATDIR = "#{ENV['PWD']}"
LIST = 'list'
INAM = 'input.dat'
CNAM = 'conf.dat'
XNAM = 'mix.dat'
MNAM = 'mie_comp.dat'
NNAM = 'nohup.out'
NRPT = 10
NPRC = 2
AMOD = [1,4,5]
NMOD = 1
refs = []
refs[0] =  "#{ENV['HOME']}/work/MODTRAN4/mie_comp/db/REFIND_9.dat #{ENV['HOME']}/work/MODTRAN4/mie_comp/db/REFIND_9.dat"
refs[1] =  "#{ENV['HOME']}/work/MODTRAN4/mie_comp/db/REFIND_4.dat #{ENV['HOME']}/work/MODTRAN4/mie_comp/db/REFIND_4.dat"
refs[2] =  "#{ENV['HOME']}/work/MODTRAN4/mie_comp/db/REFIND_3.dat #{ENV['HOME']}/work/MODTRAN4/mie_comp/db/REFIND_3.dat"

if(!File.exist?(LIST))
  $stderr.puts "No such file >>> #{LIST}"
  exit
end
dall = []
File.readlines(LIST).each{|line|
  item = line.chomp.split
  if(item.size < 1)
    $stderr.puts "Warning, item.size=#{item.size}"
    next
  end
  dall << item[0]
}

1.upto(NRPT){
  dall.each{|dnam|
    if(/^(\d\d)(\d\d)(\d\d)/ !~ dnam)
      $stderr.puts "Invalid DATE >>> #{dnam}, skip"
      next
    end
    yr = $1.to_i
    mo = $2.to_i
    dy = $3.to_i
    if(mo<1 || mo>12 || dy<1 || dy>31)
      $stderr.puts "Invalid DATE >>> yr=#{yr},mo=#{mo},dy=#{dy}, skip"
      next
    end
    caldir = "#{CALDIR}/#{dnam[0,6]}"
    if(!File.directory?(caldir))
      $stderr.puts "No such directory >>> #{caldir}, skip"
      next
    end
    iatm = 3
    if(mo>=5 && mo<=10)
      iatm = 2
    end
    ino = -1
    rno = ""
    dir = ""
    datdir = ""
    99.downto(1){|i|
      rno_tmp = sprintf("%02d",i)
      dir_tmp = "#{DATDIR}/#{dnam}/#{rno_tmp}"
      if(File.directory?(dir_tmp))
        ino = i
        rno = rno_tmp
        dir = dir_tmp
        break
      end
    }
    if(ino < 0)
      if(File.exist?("#{DATDIR}/#{dnam}/#{INAM}"))
        ino = 1
        rno = "01"
        dir = "#{DATDIR}/#{dnam}/#{rno}"
        system("mkdir -p #{dir}")
        system("cp #{DATDIR}/#{dnam}/#{INAM} #{dir}/")
        if(File.exist?("#{DATDIR}/#{dnam}/#{CNAM}"))
          system("cp #{DATDIR}/#{dnam}/#{CNAM} #{dir}/")
        end
      else
        $stderr.puts "No data for #{dnam}, skip!"
        next
      end
    end
    inam = "#{dir}/#{INAM}"
    cnam = "#{dir}/#{CNAM}"
    mnam = "#{dir}/#{MNAM}"
    nnam = "#{dir}/#{NNAM}"
    if(!File.exist?(inam))
      $stderr.puts "No such file >>> #{inam}, skip"
      next
    end
    if(!File.exist?(cnam))
      system("touch #{cnam}")
    end
    if(File.exist?(nnam))
      while(File.mtime(nnam) > Time.now-600)
        if(File.exist?(mnam))
          break
        end
        $stderr.printf("Waiting   #{mnam}   #{Time.now.strftime("%Y-%m-%d %H:%M:%S")}\n")
        sleep(300)
      end
    end
    nmod = -1
    smod = 2
    if(File.exist?(mnam))
      if(!File.exist?(nnam))
        $stderr.puts "No such file >>> #{nnam}, skip"
        next
      end
      `grep opt_amod #{nnam}`.each{|line|
        if(/opt_amod\s*:\s*(\d+)/ =~ line)
          nmod = $1.to_i
        end
      }
      if(nmod < 0)
        $stderr.puts "Warning, nmod not found >>> #{nnam}"
        nmod = NMOD
      elsif(nmod >= AMOD.size)
        $stderr.puts "Error, atmospheric model out of range >>> #{nmod}, skip (#{dnam})"
        next
      end
      lines = File.readlines(mnam)
      if(lines.size != 3)
        $stderr.puts "Error, lines.size=#{lines.size} >>> #{mnam}, skip"
        next
      end
      ino += 1
      rno = sprintf("%02d",ino)
      datdir = "#{DATDIR}/#{dnam}/#{rno}"
      if(File.exist?(datdir))
        $stderr.puts "File exist >>> #{datdir}, skip"
        next
      else
        system("mkdir -p #{datdir}")
      end
      if(!File.directory?(datdir))
        $stderr.puts "No such directory >>> #{datdir}, skip"
        next
      end
      err = 0
      File.open("#{datdir}/#{XNAM}",'w'){|fp|
        (0...lines.size).each{|i|
          item = lines[i].chomp.split
          if(item.size != 5)
            $stderr.puts "Error, item.size=#{item.size} >>> #{mnam}, skip"
            err = 1
            break
          end
          mixr = item[0]
          lmod = item[1]
          lsgm = item[2]
          fp.printf("%s %s %s %s\n",mixr,lmod,lsgm,refs[i])
        }
      }
      if(err != 0)
        next
      end
      File.open("#{datdir}/#{CNAM}",'w'){|fp|
        fcmp = 0
        File.readlines(cnam).each{|line|
          flag = 1
          item = line.chomp.split
          if(item.size < 1)
            # do nothing
          elsif(/crc_fcmp/ =~ item[0])
            fp.printf("crc_fcmp #{XNAM} 1 2 1.0e3\n")
            fcmp |= 0x01
            flag = 0
          elsif(/opt_fcmp/ =~ item[0])
            fp.printf("opt_fcmp #{XNAM} 1 2 1.0e3\n")
            fcmp |= 0x02
            flag = 0
          elsif(/prf_pres_gnd/ =~ item[0])
            smod = 1
          elsif(/prf_temp_gnd/ =~ item[0])
            smod = 1
          elsif(/prf_fhaz/ =~ item[0])
            if(item.size < 3)
              $stderr.puts "Error, item.size=#{item.size} >>> #{line}, skip"
              flag = 0
            else
              system("cp #{dir}/#{item[2]} #{datdir}/")
              smod = 1
            end
          end
          if(flag == 1)
            fp.puts line
          end
        }
        if((fcmp&0x01) < 1)
          fp.printf("crc_fcmp #{XNAM} 1 2 1.0e3\n")
        end
        if((fcmp&0x02) < 1)
          fp.printf("opt_fcmp #{XNAM} 1 2 1.0e3\n")
        end
      }
      system("cp #{inam} #{datdir}/")
    else
      nmod = NMOD
      File.readlines(cnam).each{|line|
        item = line.chomp.split
        if(item.size < 1)
          next
        elsif(/prf_pres_gnd/ =~ item[0])
          smod = 1
        elsif(/prf_temp_gnd/ =~ item[0])
          smod = 1
        elsif(/prf_fhaz/ =~ item[0])
          smod = 1
        end
      }
      datdir = dir
    end
    $stderr.flush
    flag = 0
    while(flag == 0)
      np = 0
      status = `ps -C aeros_mix -o pid,etime,cmd`
      status.each{|p|
        item = p.chomp.split
        if(item.size < 3)
          next
        end
        if(item[2] == 'aeros_mix')
          np += 1
        end
      }
      if(np < NPRC)
        system("auto_submit.sh #{datdir} #{caldir} #{smod} #{iatm} #{AMOD[nmod]} 2>/dev/null")
        $stderr.printf("Submitted #{dnam}/#{rno} @ #{Time.now.strftime("%Y-%m-%d %H:%M:%S")}\n")
        flag = 1
        sleep(1)
      else
        $stderr.printf("Waiting   #{dnam}/#{rno}   #{Time.now.strftime("%Y-%m-%d %H:%M:%S")}\n")
        sleep(300)
      end
    end
  }
}
