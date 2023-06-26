kill $(ps aux | grep 'bash scripts/open_bern2_gpu.sh' | grep -v 'grep' | awk '{print $2}');
bash BERN2/scripts/stop_bern2.sh;
kill $(fuser 8888/tcp | cut -f 2 );
