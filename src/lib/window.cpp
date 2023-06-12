/*
 *   Copyright (C) 2023 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of diffusion-problem
 *
 *   diffusion-problem is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   diffusion-problem is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with diffusion-problem.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <SFML/Config.hpp>
#include <SFML/Graphics/Texture.hpp>
#include <SFML/Graphics/Sprite.hpp>
#include "lib/window.hpp"
#include "lib/print.hpp"

namespace dplib{

void Window::save_image(){
    sf::Vector2u windowSize = window.getSize();
    sf::Texture texture;
    texture.create(windowSize.x, windowSize.y);
    texture.update(window);
    sf::Image screenshot = texture.copyToImage();
    screenshot.saveToFile(this->name + ".png");
}

Window::Window(size_t window_width, size_t window_height, size_t mesh_width, size_t mesh_height, std::string name):
    window_width(window_width), window_height(window_height), W(mesh_width), H(mesh_height),
    window(sf::VideoMode(window_width, window_height), "diffusion-problem - " + name),
    pixels(W*H*4, 255), name("diffusion-problem - " + name){

    auto resolution = sf::VideoMode::getDesktopMode();

    if(!this->font.loadFromFile("assets/LiberationSans-Regular.ttf")){
        print_line("ERROR: font loading failed.");
        exit(EXIT_FAILURE);
    }
    this->text_max.setFont(this->font);
    this->text_max.setCharacterSize(30);
    this->text_max.setFillColor(sf::Color::Black);
    this->text_max.setStyle(sf::Text::Bold);

    this->text_min.setFont(this->font);
    this->text_min.setCharacterSize(30);
    this->text_min.setFillColor(sf::Color::Black);
    this->text_min.setStyle(sf::Text::Bold);

    this->legend_width = 30;
    this->legend_height = 400;

    img.create(W, H);
    this->sprite.setTexture(this->img);
    switch(this->colormap){
        case ColorMap::GRAYSCALE:
            this->legend_grayscale();
            break;
        case ColorMap::HSV:
            // TODO
            break;
    }

    size_t offset = (window_width - W - legend_width)/3;

    //sprite.setPosition(sf::Vector2f(window_width/2-W/2, window_height/2-H/2));
    sprite.setPosition(sf::Vector2f(2*offset + legend_width, window_height/2-H/2));
    this->legend.setPosition(sf::Vector2f(offset, window_height/2-legend_height/2));
    window.setPosition(sf::Vector2i((resolution.width - window_width)/2, (resolution.height - window_height)/2));
}

void Window::update(const std::vector<double>& mesh){
    const double min_val = *std::min_element(mesh.begin(), mesh.end());
    const double max_val = *std::max_element(mesh.begin(), mesh.end());
    std::cout << min_val << " " << max_val << std::endl;
    this->update(mesh, min_val, max_val);
}
void Window::update(const std::vector<double>& mesh, const double min_x, const double max_x){
    switch(this->colormap){
        case ColorMap::GRAYSCALE:
            this->to_grayscale(mesh, min_x, max_x);
            break;
        case ColorMap::HSV:
            this->to_hsv(mesh, min_x, max_x);
            break;
    }
    this->update(min_x, max_x);
}

void Window::update(){
    sf::Event event;
    size_t offset = (window_width - W - legend_width)/3;
    while (window.pollEvent(event)){
        if (event.type == sf::Event::Closed){
            this->save_image();
            window.close();
        }
        if(event.type == sf::Event::Resized){
            sf::FloatRect view(0, 0, event.size.width, event.size.height);
            window.setView(sf::View(view));
            window_width = event.size.width;
     
            window_height = event.size.height;

            sprite.setPosition(sf::Vector2f(2*offset + legend_width, window_height/2-H/2));
            this->legend.setPosition(sf::Vector2f(offset, window_height/2-legend_height/2));
        }
    }
}

void Window::update(const double min_x, const double max_x){
    std::stringstream stream;
    stream << std::fixed << std::setprecision(3) << min_x;
    std::string s = stream.str();
    this->text_min.setString(s);
    stream.str(std::string());
    stream << std::fixed << std::setprecision(3) << max_x;
    s = stream.str();
    this->text_max.setString(s);

    size_t offset = (window_width - W - legend_width)/3;
    this->text_max.setPosition(offset + legend_width/2 - text_max.getGlobalBounds().width/2, window_height/2 - H/2 - text_max.getGlobalBounds().height - 16);
    this->text_min.setPosition(offset + legend_width/2 - text_min.getGlobalBounds().width/2, window_height/2 + H/2);
    window.clear(sf::Color(201,190,210));

    img.update(this->pixels.data());

    window.draw(sprite);
    window.draw(legend);
    window.draw(text_max);
    window.draw(text_min);
    window.display();
}

void Window::to_grayscale(const std::vector<double>& mesh, const double min_x, const double max_x){
    #pragma omp parallel for
    for(size_t i = 0; i < W*H; ++i){
        const double x = mesh[i];

        const sf::Uint8 rgb = static_cast<sf::Uint8>(255*((max_x - x)/(max_x - min_x)));
        this->pixels[i*4+0] = rgb;
        this->pixels[i*4+1] = rgb;
        this->pixels[i*4+2] = rgb;
    }
}

void Window::to_hsv(const std::vector<double>& mesh, const double min_x, const double max_x){
    const double R[]{1, 1, 0, 0, 0, 0};
    const double G[]{0, 1, 1, 1, 0, 0};
    const double B[]{0, 0, 0, 1, 1, 0};
    const size_t L = 4;
    #pragma omp parallel for
    for(size_t i = 0; i < W*H;  ++i){
        const double x = mesh[i];

        double norm = (max_x - x)/(max_x - min_x);
        const size_t low = std::floor(norm*L);
        const double rem = norm*L - low;
        const double r = R[low] + rem*(R[low+1]-R[low]);
        const double g = G[low] + rem*(G[low+1]-G[low]);
        const double b = B[low] + rem*(B[low+1]-B[low]);

        this->pixels[i*4+0] = static_cast<sf::Uint8>(255*r);
        this->pixels[i*4+1] = static_cast<sf::Uint8>(255*g);
        this->pixels[i*4+2] = static_cast<sf::Uint8>(255*b);
    }
}

void Window::legend_grayscale(){
    std::vector<sf::Uint8> px(legend_width*legend_height*4, 255);
    this->legend_img.create(legend_width, legend_height);
    #pragma omp parallel for
    for(size_t y = 0; y < legend_height; ++y){
        const sf::Uint8 p = 255.0*static_cast<double>(y)/legend_height;
        for(size_t x = 0; x < legend_width; ++x){
            for(size_t d = 0; d < 3; ++d){
                px[y*legend_width*4 + x*4 + d] = p;
            }
        }
    }
    this->legend_img.update(px.data());
    this->legend.setTexture(this->legend_img);
}

}
