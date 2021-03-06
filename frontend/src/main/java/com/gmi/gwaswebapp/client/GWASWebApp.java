package com.gmi.gwaswebapp.client;

import com.gmi.gwaswebapp.client.gin.ClientGinjector;
import com.gmi.gwaswebapp.client.resources.MyResources;
import com.google.gwt.core.client.EntryPoint;
import com.google.gwt.core.client.GWT;
import com.gwtplatform.mvp.client.DelayedBindRegistry;

public class GWASWebApp implements EntryPoint {
	
	private final ClientGinjector ginjector = GWT.create(ClientGinjector.class);

	@Override
	public void onModuleLoad() {
		DelayedBindRegistry.bind(ginjector);
		ginjector.getPlaceManager().revealCurrentPlace();
		MyResources resource = ginjector.getResource();
		resource.style().ensureInjected();
	}
}
