package statisticMerger;

import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.nio.file.StandardWatchEventKinds;
import java.nio.file.WatchEvent;
import java.nio.file.WatchKey;
import java.nio.file.WatchService;
import java.util.concurrent.TimeUnit;

public class FileWatcher {
	private boolean foundFile = false;
	
	public FileWatcher(Path monitorFile, int maxTimeoutinMs) {
		try {
		// ensure that the file is really there
		WatchService watcher = FileSystems.getDefault().newWatchService();
		Path parent = monitorFile.getParent();
		parent.register(watcher, StandardWatchEventKinds.ENTRY_MODIFY);
		int waitTime = 0;
		while(waitTime < maxTimeoutinMs && !this.foundFile) {
			waitTime = waitTime + 25;
	        WatchKey key;
	        try { key = watcher.poll(25, TimeUnit.MILLISECONDS); }
	        catch (InterruptedException e) { return; }
	        if (key == null) { continue; }

	        for(WatchEvent<?> event : key.pollEvents()) {
	            @SuppressWarnings("unchecked")
				WatchEvent<Path> ev = (WatchEvent<Path>) event;
	            Path filename = ev.context();
	            if (filename.endsWith(monitorFile.getFileName())) {
	            	this.foundFile = true;
	            	break;
	            }
	        }
	    }
		// close the watch servide
		watcher.close();
		} catch(Exception e) { e.printStackTrace(); }
	}
	
	public boolean wasSuccessfull() {
		return this.foundFile;
	}
}
