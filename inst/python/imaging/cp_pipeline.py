import matplotlib.pyplot as plt
import skimage.io
import cellprofiler_core.image
import cellprofiler_core.measurement
import cellprofiler_core.object
import cellprofiler_core.pipeline
import cellprofiler_core.preferences
import cellprofiler_core.workspace

def run_pipeline(pipeline_filename, image_dict):
    cellprofiler_core.preferences.set_headless()
    
    # Create and load the pipeline
    pipeline = cellprofiler_core.pipeline.Pipeline()
    pipeline.load(pipeline_filename)
    
    # Create the image set, and add the image data
    image_set_list = cellprofiler_core.image.ImageSetList()
    image_set = image_set_list.get_image_set(0)
    for image_name, input_pixels in image_dict.items():
        image_set.add(image_name, cellprofiler_core.image.Image(input_pixels))
        
    # Persist the object set here (for now, see workspace TODO)
    object_set = cellprofiler_core.object.ObjectSet()

    # We can only run one group -- set the group index to 1.
    measurements = cellprofiler_core.measurement.Measurements()
    measurements.group_index = 1

    # Run the modules!
    for module in pipeline.modules():
        # Yes, we really do have to create a new workspace for each module
        # because the module attribute is required. Go team.
        workspace = cellprofiler_core.workspace.Workspace(
            image_set=image_set,
            image_set_list=image_set_list,
            measurements=measurements, 
            module=module,
            object_set=object_set, 
            pipeline=pipeline
        )
        
        module.prepare_run(workspace)
        module.run(workspace)
        module.post_run(workspace)
    
    # The workspace object has access to the measurements
    # and the image set/image set list which can be used
    # to use/view/store/whatever output data.
    return workspace

def display_img (workspace):
    segmented = workspace.image_set.get_image ('CellImage')
    original = workspace.image_set.get_image ('DNA')
    skimage.io.imshow (original.pixel_data)
    skimage.io.imshow (segmented.pixel_data, alpha=0.3)
    plt.show ()
