import time
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver import ActionChains
from io import BytesIO
from base64 import b64decode
from PIL import Image
def get_url(seq,struct,colors_string):
    base_url = "http://nibiru.tbi.univie.ac.at/forna/forna.html?id=url/name&"
    colors = string_to_colors(colors_string)
    params_dict = {'sequence': seq, 'structure': struct, 'colors': colors}
    url = f'{base_url}{"&".join(f"{k}={v}"for k, v in params_dict.items())}'
    return url
def string_to_colors(colors_string):
    r,g,b,p,y = [],[],[],[],[]
    for pos,char in enumerate(colors_string):
        if(char == 'r'):
            r.append(pos+1)
        elif(char == 'g'):
            g.append(pos+1)
        elif(char == 'b'):
            b.append(pos+1)
        elif(char == 'p'):
            p.append(pos+1)
        elif(char == 'y'):
            y.append(pos+1)
    r_str = [str(x) for x in r]
    g_str = [str(x) for x in g]
    b_str = [str(x) for x in b]
    p_str = [str(x) for x in p]
    y_str = [str(x) for x in y]
    colors_dict = {
        'red': r_str, 'green': g_str, 'blue': b_str, 'pink': p_str, 'yellow': y_str
    }
    foramt_colors = ' '.join({f"{','.join(u)}:{v}" for v,u in colors_dict.items()})
    #foramt_colors = ','.join(r_str) + ':red ' + ','.join(g_str) + ':green ' + ','.join(b_str) + ':blue ' + ','.join(p_str) + ':pink'
    return foramt_colors
def get_image(url, image_path, pseudoknots_strength, driver_path):
    chrome_options = webdriver.ChromeOptions()
    chrome_options.add_argument('--no-sandbox')
    chrome_options.add_argument('--window-size=1420,1080')
    chrome_options.add_argument('--headless')
    chrome_options.add_argument('--disable-gpu')
    chrome_options.add_argument('--disable-dev-shm-usage')
    chrome_options.add_argument("--disable-extensions")
    chrome_options.add_argument("--remote-debugging-port=9222")
    chrome_options.add_argument("--shm-size=1g")
    if driver_path != 'default':
        CHROMEDRIVERPATH = driver_path
        webdriver_service = Service(CHROMEDRIVERPATH)
        with webdriver.Chrome(service=webdriver_service, options=chrome_options) as driver:
            driver.get(url)
            driver.maximize_window()
            time.sleep(2)
            toggle_menu = driver.find_elements(by=By.XPATH, value="//button[@class = 'btn btn-default dropdown-toggle']")
            if pseudoknots_strength != 'default':
                # for tight pseudoknots, to make optional
                settings_click = toggle_menu[1].click()
                pseudoknot = driver.find_element(by=By.XPATH,value="//*[contains(text(), 'Pseudoknot strength')]")
                pseudoknot_children = pseudoknot.find_elements(by=By.CSS_SELECTOR, value='*')
                slider = pseudoknot_children[1]
                move = ActionChains(driver)
                move.click_and_hold(slider).move_by_offset(10, 0).release().perform()
                #
            try:
                center_molecule = driver.find_element(by=By.XPATH, value="//button[@data-bind = 'click: centerMolecules']").click()
            except:
                pass
            toggle_menu[2].click()
            div1 = driver.find_element(by=By.XPATH, value='//div[@class = "btn-group dropup open"]')
            div1_children = div1.find_elements(by=By.CSS_SELECTOR, value='*')
            div1_children[13].click()
            get_address = driver.find_elements(by=By.XPATH, value='//a[@download= "rna.png"]')
            image_text = get_address[1].get_attribute('href')
            im = Image.open(BytesIO(b64decode(image_text.split(',')[1])))
            im.save(image_path)
        # driver.close()
    else:
        with webdriver.Chrome(options=chrome_options) as driver:
            driver.get(url)
            driver.maximize_window()
            time.sleep(2)
            toggle_menu = driver.find_elements(by=By.XPATH, value="//button[@class = 'btn btn-default dropdown-toggle']")
            if pseudoknots_strength != 'default':
                # for tight pseudoknots, to make optional
                settings_click = toggle_menu[1].click()
                pseudoknot = driver.find_element(by=By.XPATH,value="//*[contains(text(), 'Pseudoknot strength')]")
                pseudoknot_children = pseudoknot.find_elements(by=By.CSS_SELECTOR, value='*')
                slider = pseudoknot_children[1]
                move = ActionChains(driver)
                move.click_and_hold(slider).move_by_offset(10, 0).release().perform()
                #
            try:
                center_molecule = driver.find_element(by=By.XPATH, value="//button[@data-bind = 'click: centerMolecules']").click()
            except:
                pass

            toggle_menu[2].click()
            div1 = driver.find_element(by=By.XPATH, value='//div[@class = "btn-group dropup open"]')
            div1_children = div1.find_elements(by=By.CSS_SELECTOR, value='*')
            div1_children[13].click()
            get_address = driver.find_elements(by=By.XPATH, value='//a[@download= "rna.png"]')
            image_text = get_address[1].get_attribute('href')
            im = Image.open(BytesIO(b64decode(image_text.split(',')[1])))
            im.save(image_path)
        # driver.close()


def create_image(sequence: str, structure: str, colors_string: str, image_path = "RNA_image.png", pseudoknots_strength = 'default', driver_path = 'default'):
    print("Executing...")
    print("Pseudoknots strength: " + pseudoknots_strength)
    print("Driver Path: " + driver_path)
    print("Image Path: " + image_path)
    url = get_url(sequence, structure, colors_string)
    get_image(url, image_path, pseudoknots_strength, driver_path)
    print("Done!")

if __name__ == '__main__':
    chrome_options = webdriver.ChromeOptions()
    chrome_options.add_argument('--no-sandbox')
    chrome_options.add_argument('--window-size=1420,1080')
    chrome_options.add_argument('--headless')
    chrome_options.add_argument('--disable-gpu')
    chrome_options.add_argument('--disable-dev-shm-usage')
    chrome_options.add_argument("--disable-extensions")
    chrome_options.add_argument("--remote-debugging-port=9222")
    chrome_options.add_argument("--shm-size=1g")
    driver = webdriver.Chrome(options=chrome_options)
    driver.close()
    create_image("AAAACCCUUAUUACGAUUAGUAUUGAUAUUCUAUAAAAAAUUCGCCACCAUGGGAUAUUGAUACUGAAAACCUGGCGGCAGCGCAAAAG", "...................((((((((((((((.....[[[[[....)))))))))))))).....]]]]]..................", 'rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrggggggggggggggggggggyyyyyyyyyyyyyyyyyyyyyyyyyyyyy', r"C:\Users\efimo\Desktop\IGEM 2022\COLOR RNA\finally22.png")
