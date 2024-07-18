import json
from importlib.resources import path
from pathlib import Path
from typing import List, Optional

import dataclasses
import uuid
import dataclasses
import pathlib
import tempfile
from Color_RNA.src.Color_RNA import Color_RNA
from fpdf import FPDF
from typing import List
import uuid
from PIL import Image


@dataclasses.dataclass
class SwitchResults:
    switch: str
    switch_structure: str
    switch_colors: str

    metric_score: float
    is_kissing_hairpin: bool
    regressor_on: float
    regressor_off: float
    found_epitops: bool = False


@dataclasses.dataclass
class TriggerResult:
    trigger: str
    add_what_ever_you_want: bool



@dataclasses.dataclass
class PipelineResults:
    job_id: str
    trigger_result: TriggerResult
    switch_results: List[SwitchResults]
    optimized_gene: Optional[str]



class PDF(FPDF):
     def charts(self):
        self.set_xy(40.0, 25.0)
        self.image(link='', type='', w=700/5, h=450/5)




#Efi's added function for loading json to the right format
def load_pipline_results(path):
    with open(path, 'r') as f:
        jsondict = json.loads(f.read())
        switch_results_list = []
    for d in jsondict['switch_results']:
        switch_results_list.append(SwitchResults(
            switch = d['switch'],
            switch_structure = d['switch_structure'],
            switch_colors = d['switch_colors'],
            metric_score = d['metric_score'],
            is_kissing_hairpin = d['is_kissing_hairpin'],
            regressor_on = d['regressor_on'],
            regressor_off = d['regressor_off']
        ))
    trig_res = TriggerResult(trigger = jsondict['trigger_result'],add_what_ever_you_want=0)
    pipeline_res = PipelineResults(
        job_id = jsondict['job_id'],
        trigger_result = trig_res,
        switch_results = switch_results_list,
        optimized_gene = None
    )

    return pipeline_res




def generate_pdf(pipeline_results: PipelineResults, dst_path: Path, order_by: str):
    pdf = PDF(format='A4')

    with tempfile.TemporaryDirectory() as temp_dir:
        output_images_temp_folder = pathlib.Path(temp_dir)
        pdf.add_page()

        if order_by == 'Dynamic Range':
            sorted_switch_results = sorted(pipeline_results.switch_results, key=lambda switch: switch.metric_score)
        elif order_by == 'Off':
            sorted_switch_results = sorted(pipeline_results.switch_results, key=lambda switch: switch.regressor_off)
        else:
            sorted_switch_results = sorted(pipeline_results.switch_results, key=lambda switch: switch.regressor_on)

        for idx, switch_result in enumerate(sorted_switch_results):
            switch_image_name = str(output_images_temp_folder / f"{uuid.uuid4()}.png")
            switch_image_name_jpg = str(output_images_temp_folder / f"{uuid.uuid4()}.jpg")
            Color_RNA.create_image(switch_result.switch, switch_result.switch_structure, switch_result.switch_colors, switch_image_name)
            im1 = Image.open(switch_image_name)
            rgb_im = im1.convert('RGB')
            rgb_im.save(switch_image_name_jpg)

            pos = 20 + (50 * idx)
            if idx == 0:
                pdf.set_font('Arial', 'BU', 14)
                pdf.cell(0, 0, "Results from TrigGate's Switch Generator", 0, 2, 'C')


            pdf.set_xy(10, pos-5)
            pdf.set_font('Arial', 'BU', 10)
            pdf.write(5, "Result #" +str(idx+1))

            pdf.set_xy(10, pos)
            pdf.set_font('Arial', 'B', 7)
            pdf.write(5, f"{switch_result.switch}")

            pdf.set_xy(10, pos + 15)
            pdf.set_font('Arial', size=11)
            pdf.write(5, f"Metric_score: {switch_result.metric_score}")

            pdf.ln()
            pdf.write(5, f"Is kissing hairpin: {switch_result.is_kissing_hairpin}")

            pdf.ln()
            pdf.write(5, f"Regressor on rank: {switch_result.regressor_on}")

            pdf.ln()
            pdf.write(5, f"Regressor off rank: {switch_result.regressor_off}")

            pdf.set_xy(135.0, pos + 5)
            pdf.image(switch_image_name_jpg,  link='', type='', w=700/15, h=450/15)
            pdf.line(5,pos+40,205,pos+40)
        if pipeline_results.optimized_gene is not None:
            pdf.add_page()
            pdf.set_font('Arial', 'BU', 14)
            pdf.cell(0, 0, "Translation Optimization Results", 0, 2, 'C')

            pdf.set_xy(10, 20)
            pdf.set_font('Arial',size=10)
            pdf.write(5, "Gene Optimization Results: " + pipeline_results.optimized_gene)


        pdf.output(str(pathlib.Path(dst_path) / f'{pipeline_results.job_id}.pdf'), 'F')