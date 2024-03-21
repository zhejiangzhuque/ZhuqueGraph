from PIL import Image, ImageDraw,ImageFont
import json
import os

class Wall:
    def __init__(self, images, width, height, coordinates) -> None:
        self.images = images
        self.width = width
        self.height = height
        self.coordinates = coordinates

    def draw(self):
        for coordinate in self.coordinates:
            x, y = coordinate
            for image in self.images:
                image.putpixel((x, y), (127, 127, 127, 1))  # 灰色

class View:
    def __init__(self, config_path, info_path, save_path) -> None:
        # 图像大小
        self.background_color = (255, 255, 255)
        self.scale_coef = 10
        self.save_path = save_path
        self.config_path = config_path
        self.info_path = info_path
        with open(os.path.join(config_path), "r") as f:
            config = json.load(f)
            self.width = config["width"]
            self.height = config["height"]
        self.coordinates = []
        self.wall_coordinates = []
        self.wall_flag = False
        self.agent_flag = False
        self.num_coordinates = []
        self.agent_counts = []

    def get_coordinates(self):
        # 解析信息行和图形坐标
        with open(os.path.join(self.info_path), "r") as file:
            lines = file.readlines()
            for line in lines:
                elements = line.strip().split(" ")
                if elements[0] == "F" or elements[0] == "W":
                    if elements[0] == "F":
                        num_coordinate = int(elements[1])
                        self.num_coordinates.append(num_coordinate)
                        self.agent_counts.append([0,0])
                        self.agent_flag = True
                        self.wall_flag = False
                        agent_count = 0
                    elif elements[0] == "W":
                        self.agent_flag = False
                        self.wall_flag = True
                else:
                    if elements[0] == "0" and len(elements) == 4:
                        continue
                    if self.agent_flag:
                        color = (255, 0, 0) if int(elements[0]) < 64 else (0, 0, 255)
                        angle = int(elements[2])
                        x = int(elements[3])
                        y = int(elements[4])
                        self.coordinates.append((x, y, color))
                        if int(elements[0]) < 64:
                            self.agent_counts[-1][0] += 1
                        else:
                            self.agent_counts[-1][1] += 1
                    if self.wall_flag:
                        x = int(elements[0])
                        y = int(elements[1])
                        self.wall_coordinates.append((x, y))
            return agent_count

    def draw(self):
        # 创建图像对象
        agent_count = self.get_coordinates()
        images = [Image.new("RGB", (self.width, self.height), self.background_color) for _ in range(len(self.num_coordinates))]
        wall = Wall(images, self.width, self.height, self.wall_coordinates)
        wall.draw()

        # 绘制坐标点
        idx = 0
        for image_idx, num_coordinate in enumerate(self.num_coordinates):
            for coordinate_idx in range(idx, num_coordinate + idx):
                x, y, color = self.coordinates[coordinate_idx]
                images[image_idx].putpixel((x, y), color)

            idx += num_coordinate

        # 在图像上显示智能体数量
        for i, image in enumerate(images):
            draw = ImageDraw.Draw(image)

        # 图形放大
        resized_images = [image.resize((self.width*self.scale_coef, self.height*self.scale_coef), resample=Image.NEAREST) for image in images]
        # 保存为GIF动画
        for i, image in enumerate(resized_images):
            draw = ImageDraw.Draw(image)
            count_text = f"Red Num:  {self.agent_counts[i][0]}\nBlue Num: {self.agent_counts[i][1]}"
            draw.text((20, 20), count_text, fill=(0, 0, 0))
        resized_images[0].save(self.save_path, save_all=True, append_images=resized_images[1:], optimize=False, duration=200, loop=0)

if __name__ == '__main__':
    config_path = r"/dev/workspace/LSC/rl/algorithm/config.json"
    info_path = r"/dev/workspace/LSC/rl/algorithm/coordinates.txt"
    save_path = r"/dev/workspace/LSC/rl/algorithm/coordinates.gif"
    view = View(config_path, info_path, save_path)
    view.draw()
