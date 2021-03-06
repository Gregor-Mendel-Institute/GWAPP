package com.gmi.gwaswebapp.client.mvp.analysis;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


import com.gmi.gwaswebapp.client.NameTokens;
import com.gmi.gwaswebapp.client.dto.Analysis;
import com.gmi.gwaswebapp.client.dto.BaseModel;
import com.gmi.gwaswebapp.client.dto.BaseModelKeyProvider;
import com.gmi.gwaswebapp.client.dto.Dataset;
import com.gmi.gwaswebapp.client.dto.Phenotype;
import com.gmi.gwaswebapp.client.dto.Transformation;
import com.gmi.gwaswebapp.client.mvp.analysis.AnalysisPresenter.MyView;
import com.gmi.gwaswebapp.client.ui.SlidingPanel;
import com.gmi.gwaswebapp.client.ui.SlidingPanel.DIRECTION;
import com.google.gwt.core.client.GWT;
import com.google.gwt.event.logical.shared.SelectionEvent;
import com.google.gwt.event.logical.shared.SelectionHandler;
import com.google.gwt.uibinder.client.UiBinder;
import com.google.gwt.user.client.ui.HTMLPanel;
import com.google.gwt.user.client.ui.Hyperlink;
import com.google.gwt.user.client.ui.SourcesTreeEvents;
import com.google.gwt.user.client.ui.Tree;
import com.google.gwt.user.client.ui.TreeItem;
import com.google.gwt.user.client.ui.Widget;
import com.google.inject.Inject;
import com.gwtplatform.mvp.client.ViewWithUiHandlers;
import com.gwtplatform.mvp.client.proxy.PlaceManager;
import com.gwtplatform.mvp.client.proxy.PlaceRequest;
import com.google.gwt.uibinder.client.UiField;
import com.google.gwt.user.cellview.client.CellTree;
import com.google.gwt.user.cellview.client.TreeNode;
import com.google.gwt.view.client.SelectionChangeEvent;
import com.google.gwt.view.client.SingleSelectionModel;


public class AnalysisView extends ViewWithUiHandlers<AnalysisUiHandlers> implements MyView {

	private static AnalysisViewUiBinder uiBinder = GWT
			.create(AnalysisViewUiBinder.class);
//	@UiField(provided=true) CellTree resultsTree;
	@UiField Tree tree;

	interface AnalysisViewUiBinder extends UiBinder<Widget, AnalysisView> {
	}
	
	@UiField SlidingPanel contentContainer;
	//@UiField HTMLPanel contentContainer2;
	
	private Widget widget;
	private SingleSelectionModel<BaseModel> selectionModel = null;
	private AnalysisTreeViewModel analysisTreeViewModel = null;
	private final PlaceManager placeManager;

	@Inject
	public AnalysisView(final PlaceManager placeManager) {
		this.placeManager = placeManager;
		/*CellTree.Resources res = GWT.create(CellTree.BasicResources.class);
		selectionModel  = new SingleSelectionModel<BaseModel>(BaseModel.KEY_PROVIDER);
		selectionModel.addSelectionChangeHandler(new SelectionChangeEvent.Handler() {
			@Override
			public void onSelectionChange(SelectionChangeEvent event) {
				Object object = selectionModel.getSelectedObject();
				if (object instanceof Analysis)
				{
					Analysis analysis = (Analysis)selectionModel.getSelectedObject();
					if (getUiHandlers() != null && analysis != null)
						getUiHandlers().onLoadGWAS(analysis);
				}
				else if (object instanceof Dataset) {
					Dataset dataset = (Dataset)selectionModel.getSelectedObject();
					if (getUiHandlers() != null && dataset != null) 
						getUiHandlers().onLoadDataset(dataset);
				}
				else if (object instanceof Transformation)
				{
					Transformation transformation  = (Transformation)selectionModel.getSelectedObject();
					if (getUiHandlers() != null && transformation != null)
						getUiHandlers().onLoadTransformation(transformation);
				}
				else if (object instanceof Phenotype)
				{
					Phenotype phenotype = (Phenotype)selectionModel.getSelectedObject();
					if (getUiHandlers() != null)
						getUiHandlers().onLoadPhenotype(phenotype);
				}
			}
		});
		analysisTreeViewModel = new AnalysisTreeViewModel(new ArrayList<Phenotype>(),selectionModel);
		resultsTree = new CellTree(analysisTreeViewModel,null,res);
		resultsTree.setAnimationEnabled(true);*/
		widget = uiBinder.createAndBindUi(this);
		contentContainer.setDirection(DIRECTION.VERTICAL);
		tree.addSelectionHandler(new SelectionHandler<TreeItem>() {
			
			@Override
			public void onSelection(SelectionEvent<TreeItem> event) {
				Object object = event.getSelectedItem().getUserObject();
				if (object instanceof Analysis)
				{
					Analysis analysis = (Analysis)object;
					if (getUiHandlers() != null && analysis != null)
						getUiHandlers().onLoadGWAS(analysis);
				}
				else if (object instanceof Dataset) {
					Dataset dataset = (Dataset)object;
					if (getUiHandlers() != null && dataset != null) 
						getUiHandlers().onLoadDataset(dataset);
				}
				else if (object instanceof Transformation)
				{
					Transformation transformation  = (Transformation)object;
					if (getUiHandlers() != null && transformation != null)
						getUiHandlers().onLoadTransformation(transformation);
				}
				else if (object instanceof Phenotype)
				{
					Phenotype phenotype = (Phenotype)object;
					if (getUiHandlers() != null)
						getUiHandlers().onLoadPhenotype(phenotype);
				}
				
			}
		});
	}

	@Override
	public Widget asWidget() {
		return widget;
	}

	public void initTreeItems(List<Phenotype> phenotypes) {
		
		//analysisTreeViewModel.updateTree(phenotypes);
		tree.clear();
		for (Phenotype phenotype: phenotypes) {
			TreeItem phenotypeItem = new TreeItem(phenotype.getName());
			phenotypeItem.setUserObject(phenotype);
			for (Dataset dataset: phenotype.getDatasets()) {
				TreeItem datasetItem = new TreeItem(dataset.getName());
				datasetItem.setUserObject(dataset);
				for (Transformation transformation:dataset.getTransformations()) {
					TreeItem transformationItem = new TreeItem(transformation.getName());
					transformationItem.setUserObject(transformation);
					for (Analysis analysis: transformation.getAnalysisMethods()) {
						TreeItem analysisItem= new TreeItem(analysis.getName());
						analysisItem.setUserObject(analysis);
						transformationItem.addItem(analysisItem);
					}
					datasetItem.addItem(transformationItem);
				}
				phenotypeItem.addItem(datasetItem);
			}
			tree.addItem(phenotypeItem);
		}
	}
	
	@Override
	public void setInSlot(Object slot,Widget content) {
		if (slot == AnalysisPresenter.TYPE_SetPhenotypeListContent) {
			setPhenotypeListContent(content);
		}	
		else if (slot==AnalysisPresenter.TYPE_SetPhenotypeDetailContent) {
			setPhenotypeDetailContent(content);
		}
		else if (slot == AnalysisPresenter.TYPE_SetResultDetailContent)
			setResultDetailContent(content);
		else
			super.setInSlot(slot,content);
	}
	
	protected void setPhenotypeListContent(Widget content) {
		if (content != null)
			contentContainer.setWidget(content);
	}
	
	protected void setPhenotypeDetailContent(Widget content) {
		if (content != null)
			contentContainer.setWidget(content);
	}
	
	protected void setResultDetailContent(Widget content) {
		if (content != null)
			contentContainer.setWidget(content);
	}

	@Override
	public SingleSelectionModel<BaseModel> getTreeSelectionModel() {
		return selectionModel;
	}

	
	@Override
	public void setSelectedTreeItem(BaseModel itemToSelect) {
		if (itemToSelect == null) {
			tree.setSelectedItem(null);
			return;
		}
		Iterator<TreeItem> iterator = null;
		BaseModel item = null;
		if (tree.getSelectedItem() != null) {
			item = (BaseModel)tree.getSelectedItem().getUserObject();
			if (BaseModel.KEY_PROVIDER.getKey(item).equals(BaseModel.KEY_PROVIDER.getKey(itemToSelect))) {
				tree.ensureSelectedItemVisible();
				return;
			}
				
		}
		iterator = tree.treeItemIterator();
		while(iterator.hasNext()) {
			TreeItem treeItem = iterator.next();
			item = (BaseModel)treeItem.getUserObject();
			if (BaseModel.KEY_PROVIDER.getKey(item).equals(BaseModel.KEY_PROVIDER.getKey(itemToSelect))) {
				tree.setSelectedItem(treeItem,false);
				tree.ensureSelectedItemVisible();
				return;
			}
		}
	}
	
	@Override
	public boolean expandTree(BaseModel selectedItem,TreeNode node) {
		
		/*if (node == null)
			node = resultsTree.getRootTreeNode();
		int childCount = node.getChildCount();
		if (node.getValue() == selectedItem) {
			return true;
		}
		for (int i=0;i<childCount;i++) {
			if (node.isChildLeaf(i)) 
			{
				if (node.getChildValue(i) == selectedItem)
					return true;
				else
					node.setChildOpen(i, false);
			}
			else if (!expandTree(selectedItem, node.setChildOpen(i, true)))
				node.setChildOpen(i, false);
			else 
				return true;
		}*/
		return false;
	}

}
